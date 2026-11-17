"""
2D sagittal Anterior Chamber Boussinesq thermal-buoyancy fluid dynamics and
intracameral pigment / drug advection-diffusion PDE solver.

Models the physiological anterior chamber convection current driven by the
temperature differential between the warm iris/pupil surface (~37.0 C) and
the cooler corneal endothelium (~34.0 C) exposed to ambient air, using the
2D steady/pseudo-transient streamfunction-vorticity formulation:
    nabla^2 psi = -omega
    u = d(psi)/dy,   v = -d(psi)/dx
    d(omega)/dt + u * d(omega)/dx + v * d(omega)/dy = nu * nabla^2 omega + g * beta * dT/dx
    dT/dt + u * dT/dx + v * dT/dy = alpha_th * nabla^2 T
    dc/dt + u * dc/dx + (v - v_settle) * dc/dy = D_c * nabla^2 c

where `c(x, y)` tracks intracameral pigment granules (Krukenberg spindle
deposition on the inferior corneal endothelium) or topically absorbed drug.
"""

"""
    solve_anterior_convection_2d(anatomy; nx=33, ny=49, t_cornea_C=34.0, t_iris_C=37.0, n_time_steps=180, dt_s=0.05, pigment_settling_mm_s=0.025, iop_mmHg=15.0)

Solve the coupled 2D sagittal Boussinesq Navier-Stokes (streamfunction-vorticity),
energy (temperature), and scalar tracer advection-diffusion equations inside the
anterior chamber cross-section (`x` in `[0, depth]`, `y` in `[0, height]`).

Returns an `AnteriorChamberPDEState` containing velocity, streamfunction,
temperature, tracer deposition fields, Nusselt number, and convergence metrics.
"""
function solve_anterior_convection_2d(
    anatomy::OcularAnatomy;
    nx::Int = 33,
    ny::Int = 49,
    t_cornea_C::Float64 = 34.0,
    t_iris_C::Float64 = 37.0,
    n_time_steps::Int = 180,
    dt_s::Float64 = 0.05,
    pigment_settling_mm_s::Float64 = 0.025,
    iop_mmHg::Float64 = 15.0,
)
    nx >= 11 && ny >= 15 || throw(ArgumentError("Grid dimensions must be at least 11x15"))
    W = anatomy.anterior_chamber_depth_mm
    H = anatomy.anterior_chamber_height_mm
    dx = W / (nx - 1)
    dy = H / (ny - 1)

    x_mm = collect(range(0.0, W; length = nx))
    y_mm = collect(range(0.0, H; length = ny))

    # Aqueous humor physical properties in mm-s-C units at ~35.5 C
    nu = 0.72          # kinematic viscosity (mm^2/s)
    alpha_th = 0.148   # thermal diffusivity (mm^2/s)
    D_c = 0.012        # effective tracer/pigment diffusivity (mm^2/s)
    g_beta = 9810.0 * 3.4e-4 # g * thermal expansion (mm / (s^2 * K))

    # Slight dampening of convection under severe ocular hypertension / shallowing
    visc_mod = 1.0 + 0.008 * max(iop_mmHg - 15.0, 0.0)
    nu_eff = nu * visc_mod

    T = zeros(Float64, nx, ny)
    T_new = zeros(Float64, nx, ny)
    psi = zeros(Float64, nx, ny)
    omega = zeros(Float64, nx, ny)
    omega_new = zeros(Float64, nx, ny)
    u = zeros(Float64, nx, ny)
    v = zeros(Float64, nx, ny)
    C = zeros(Float64, nx, ny)
    C_new = zeros(Float64, nx, ny)

    # Linear initial temperature profile: x=0 cornea (cool), x=W iris (warm)
    for j in 1:ny, i in 1:nx
        frac = (i - 1) / (nx - 1)
        T[i, j] = t_cornea_C + frac * (t_iris_C - t_cornea_C)
    end
    # Initial pigment/drug distribution in central pupillary/anterior chamber core
    for j in 1:ny, i in 1:nx
        rx = (x_mm[i] - 0.50 * W) / (0.28 * W)
        ry = (y_mm[j] - 0.50 * H) / (0.25 * H)
        C[i, j] = exp(-(rx^2 + ry^2))
    end
    c_init_sum = sum(C)
    if c_init_sum > 0.0
        C ./= (c_init_sum * dx * dy)
    end

    # Restrict time step for explicit diffusion/advection stability
    dt_diff = 0.20 * min(dx^2, dy^2) / max(nu_eff, alpha_th, D_c)
    dt = min(dt_s, dt_diff)

    idx2 = 1.0 / (dx * dx)
    idy2 = 1.0 / (dy * dy)
    rhs_psi = zeros(Float64, nx, ny)
    sor_iters = 0
    sor_res = 0.0
    max_res = 0.0

    for _ in 1:n_time_steps
        # 1. Solve nabla^2 psi = -omega via SOR
        rhs_psi .= .-omega
        sor_iters, sor_res = solve_poisson_sor!(psi, rhs_psi, dx, dy; omega = 1.72, tol = 1e-5, max_iter = 400)

        # 2. Recover velocity u = dpsi/dy, v = -dpsi/dx
        @inbounds for j in 2:(ny - 1), i in 2:(nx - 1)
            u[i, j] = (psi[i, j + 1] - psi[i, j - 1]) / (2.0 * dy)
            v[i, j] = -(psi[i + 1, j] - psi[i - 1, j]) / (2.0 * dx)
        end

        # 3. No-slip wall vorticity boundary conditions (Thom's formula)
        @inbounds for j in 2:(ny - 1)
            omega[1, j] = -2.0 * psi[2, j] * idx2
            omega[nx, j] = -2.0 * psi[nx - 1, j] * idx2
        end
        @inbounds for i in 2:(nx - 1)
            omega[i, 1] = -2.0 * psi[i, 2] * idy2
            omega[i, ny] = -2.0 * psi[i, ny - 1] * idy2
        end

        # 4. Advance vorticity with upwind advection + diffusion + buoyancy
        @inbounds for j in 2:(ny - 1), i in 2:(nx - 1)
            uij = u[i, j]
            vij = v[i, j]
            dw_dx = uij >= 0.0 ? (omega[i, j] - omega[i - 1, j]) / dx : (omega[i + 1, j] - omega[i, j]) / dx
            dw_dy = vij >= 0.0 ? (omega[i, j] - omega[i, j - 1]) / dy : (omega[i, j + 1] - omega[i, j]) / dy
            lap_w = (omega[i + 1, j] - 2.0 * omega[i, j] + omega[i - 1, j]) * idx2 +
                    (omega[i, j + 1] - 2.0 * omega[i, j] + omega[i, j - 1]) * idy2
            dT_dx = (T[i + 1, j] - T[i - 1, j]) / (2.0 * dx)
            omega_new[i, j] = omega[i, j] + dt * (-uij * dw_dx - vij * dw_dy + nu_eff * lap_w + g_beta * dT_dx)
        end
        omega[2:(nx - 1), 2:(ny - 1)] .= omega_new[2:(nx - 1), 2:(ny - 1)]

        # 5. Advance Temperature T with Dirichlet x-walls and adiabatic top/bottom
        @inbounds for j in 2:(ny - 1), i in 2:(nx - 1)
            uij = u[i, j]
            vij = v[i, j]
            dT_dx = uij >= 0.0 ? (T[i, j] - T[i - 1, j]) / dx : (T[i + 1, j] - T[i, j]) / dx
            dT_dy = vij >= 0.0 ? (T[i, j] - T[i, j - 1]) / dy : (T[i, j + 1] - T[i, j]) / dy
            lap_T = (T[i + 1, j] - 2.0 * T[i, j] + T[i - 1, j]) * idx2 +
                    (T[i, j + 1] - 2.0 * T[i, j] + T[i, j - 1]) * idy2
            T_new[i, j] = T[i, j] + dt * (-uij * dT_dx - vij * dT_dy + alpha_th * lap_T)
        end
        T_new[1, :] .= t_cornea_C
        T_new[nx, :] .= t_iris_C
        T_new[:, 1] .= T_new[:, 2]
        T_new[:, ny] .= T_new[:, ny - 1]
        max_res = maximum(abs.(T_new .- T)) / dt
        T .= T_new

        # 6. Conservative Finite-Volume upwind advection-diffusion for tracer C
        #    Zero-flux Neumann boundaries preserve total intracameral tracer mass!
        @inbounds for j in 2:(ny - 1), i in 2:(nx - 1)
            u_e = 0.5 * (u[i, j] + u[i + 1, j])
            u_w = 0.5 * (u[i - 1, j] + u[i, j])
            v_n = 0.5 * (v[i, j] + v[i, j + 1]) - pigment_settling_mm_s
            v_s = 0.5 * (v[i, j - 1] + v[i, j]) - pigment_settling_mm_s

            flux_e = (u_e >= 0.0 ? u_e * C[i, j] : u_e * C[i + 1, j]) - D_c * (C[i + 1, j] - C[i, j]) / dx
            flux_w = (u_w >= 0.0 ? u_w * C[i - 1, j] : u_w * C[i, j]) - D_c * (C[i, j] - C[i - 1, j]) / dx
            flux_n = (v_n >= 0.0 ? v_n * C[i, j] : v_n * C[i, j + 1]) - D_c * (C[i, j + 1] - C[i, j]) / dy
            flux_s = (v_s >= 0.0 ? v_s * C[i, j - 1] : v_s * C[i, j]) - D_c * (C[i, j] - C[i, j - 1]) / dy

            C_new[i, j] = max(0.0, C[i, j] - dt * ((flux_e - flux_w) / dx + (flux_n - flux_s) / dy))
        end
        # Enforce zero normal flux on boundaries
        C_new[1, :] .= C_new[2, :]
        C_new[nx, :] .= C_new[nx - 1, :]
        C_new[:, 1] .= C_new[:, 2]
        C_new[:, ny] .= C_new[:, ny - 1]
        C .= C_new
    end

    speed = sqrt.(u .^ 2 .+ v .^ 2)
    max_vel = maximum(speed)

    # Compute mean Nusselt number at corneal endothelium (x = 0)
    delta_T = max(abs(t_iris_C - t_cornea_C), 1e-6)
    nu_sum = 0.0
    for j in 2:(ny - 1)
        dT_dx_wall = (T[2, j] - T[1, j]) / dx
        nu_sum += (W / delta_T) * dT_dx_wall
    end
    mean_nusselt = nu_sum / (ny - 2)

    # Inferior half vs superior half tracer deposition (Krukenberg spindle metric)
    mid_j = div(ny, 2)
    inf_mass = sum(C[:, 1:mid_j])
    tot_mass = max(sum(C), 1e-12)
    inferior_ratio = inf_mass / tot_mass

    return AnteriorChamberPDEState(
        nx,
        ny,
        W,
        H,
        x_mm,
        y_mm,
        T,
        psi,
        omega,
        u,
        v,
        C,
        max_vel,
        mean_nusselt,
        inferior_ratio,
        n_time_steps,
        max_res + sor_res,
    )
end
