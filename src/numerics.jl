"""
Pure-Julia numerical integrators, elliptic PDE solvers, root finders,
nonlinear least-squares calibration, and global sensitivity analysis routines.

Implements:
- Adaptive Runge-Kutta-Fehlberg 4(5) (`rkf45_integrate`) with Cash-Karp tableau
  and PI step-size control.
- Stiff Implicit Backward Euler (`implicit_euler_integrate`) and Crank-Nicolson
  (`crank_nicolson_integrate`) with finite-difference Jacobian and damped Newton.
- Successive Over-Relaxation (SOR) elliptic Poisson solver (`solve_poisson_sor!`).
- Brent-Dekker hybrid scalar root finder (`brent_root`).
- Levenberg-Marquardt nonlinear least-squares solver (`levenberg_marquardt`).
- Welford online statistics (`WelfordAccumulator`) and Morris elementary effects
  global sensitivity analysis (`morris_sensitivity`).
"""

using LinearAlgebra
using Random
using Statistics

# ---------------------------------------------------------------------------
# Cash-Karp Runge-Kutta-Fehlberg 4(5) Tableau Coefficients
# ---------------------------------------------------------------------------
const _CK_A = (0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0)

const _CK_B21 = 1.0 / 5.0
const _CK_B31 = 3.0 / 40.0
const _CK_B32 = 9.0 / 40.0
const _CK_B41 = 3.0 / 10.0
const _CK_B42 = -9.0 / 10.0
const _CK_B43 = 6.0 / 5.0
const _CK_B51 = -11.0 / 54.0
const _CK_B52 = 5.0 / 2.0
const _CK_B53 = -70.0 / 27.0
const _CK_B54 = 35.0 / 27.0
const _CK_B61 = 1631.0 / 55296.0
const _CK_B62 = 175.0 / 512.0
const _CK_B63 = 575.0 / 13824.0
const _CK_B64 = 44275.0 / 110592.0
const _CK_B65 = 253.0 / 4096.0

# 5th-order weights
const _CK_C1 = 37.0 / 378.0
const _CK_C3 = 250.0 / 621.0
const _CK_C4 = 125.0 / 594.0
const _CK_C6 = 512.0 / 1771.0

# Error difference (5th - 4th order)
const _CK_DC1 = _CK_C1 - 2825.0 / 27648.0
const _CK_DC3 = _CK_C3 - 18575.0 / 48384.0
const _CK_DC4 = _CK_C4 - 13525.0 / 55296.0
const _CK_DC5 = -277.0 / 14336.0
const _CK_DC6 = _CK_C6 - 1.0 / 4.0

"""
    rkf45_integrate(rhs!, y0, tspan; rtol=1e-6, atol=1e-8, h0=0.0, hmin=1e-9, hmax=1.0, max_steps=200_000, sample_times=nothing)

Integrate `dy/dt = rhs!(dy, y, t)` over `tspan = (t0, tf)` using the adaptive
5th-order Cash-Karp Runge-Kutta-Fehlberg embedded pair.

If `sample_times` is provided, values are linearly/hermite-interpolated onto those
exact timestamps; otherwise all accepted adaptive steps are returned.
Returns `(times::Vector{Float64}, states::Matrix{Float64})` where `states[:, k]`
is the state vector at `times[k]`.
"""
function rkf45_integrate(
    rhs!::Function,
    y0::AbstractVector{Float64},
    tspan::Tuple{Float64,Float64};
    rtol::Float64 = 1e-6,
    atol::Float64 = 1e-8,
    h0::Float64 = 0.0,
    hmin::Float64 = 1e-9,
    hmax::Float64 = 1.0,
    max_steps::Int = 200_000,
    sample_times::Union{Nothing,AbstractVector{Float64}} = nothing,
)
    t0, tf = Float64(tspan[1]), Float64(tspan[2])
    tf >= t0 || throw(ArgumentError("tspan must satisfy tf >= t0"))
    n = length(y0)
    y = copy(Float64.(y0))

    if tf == t0
        return [t0], reshape(copy(y), n, 1)
    end

    h = h0 > 0.0 ? min(h0, tf - t0, hmax) : min((tf - t0) / 50.0, hmax, 0.05)
    h = max(h, hmin)

    k1 = zeros(Float64, n)
    k2 = zeros(Float64, n)
    k3 = zeros(Float64, n)
    k4 = zeros(Float64, n)
    k5 = zeros(Float64, n)
    k6 = zeros(Float64, n)
    ytmp = zeros(Float64, n)
    ynew = zeros(Float64, n)
    yerr = zeros(Float64, n)

    raw_t = Float64[t0]
    raw_y = Vector{Float64}[copy(y)]
    raw_dy = Vector{Float64}[]
    dy0 = zeros(Float64, n)
    rhs!(dy0, y, t0)
    push!(raw_dy, copy(dy0))

    t = t0
    steps = 0
    while t < tf - 1e-13
        steps += 1
        steps <= max_steps || error("rkf45_integrate exceeded max_steps=$max_steps at t=$t")

        if t + h > tf
            h = tf - t
        end

        rhs!(k1, y, t)
        @inbounds for i in 1:n
            ytmp[i] = y[i] + h * _CK_B21 * k1[i]
        end
        rhs!(k2, ytmp, t + _CK_A[2] * h)

        @inbounds for i in 1:n
            ytmp[i] = y[i] + h * (_CK_B31 * k1[i] + _CK_B32 * k2[i])
        end
        rhs!(k3, ytmp, t + _CK_A[3] * h)

        @inbounds for i in 1:n
            ytmp[i] = y[i] + h * (_CK_B41 * k1[i] + _CK_B42 * k2[i] + _CK_B43 * k3[i])
        end
        rhs!(k4, ytmp, t + _CK_A[4] * h)

        @inbounds for i in 1:n
            ytmp[i] = y[i] + h * (_CK_B51 * k1[i] + _CK_B52 * k2[i] + _CK_B53 * k3[i] + _CK_B54 * k4[i])
        end
        rhs!(k5, ytmp, t + _CK_A[5] * h)

        @inbounds for i in 1:n
            ytmp[i] = y[i] + h * (
                _CK_B61 * k1[i] + _CK_B62 * k2[i] + _CK_B63 * k3[i] + _CK_B64 * k4[i] + _CK_B65 * k5[i]
            )
        end
        rhs!(k6, ytmp, t + _CK_A[6] * h)

        err_norm = 0.0
        @inbounds for i in 1:n
            ynew[i] = y[i] + h * (_CK_C1 * k1[i] + _CK_C3 * k3[i] + _CK_C4 * k4[i] + _CK_C6 * k6[i])
            yerr[i] = h * (_CK_DC1 * k1[i] + _CK_DC3 * k3[i] + _CK_DC4 * k4[i] + _CK_DC5 * k5[i] + _CK_DC6 * k6[i])
            sc = atol + rtol * max(abs(y[i]), abs(ynew[i]))
            r = abs(yerr[i]) / sc
            err_norm = max(err_norm, r)
        end

        if err_norm <= 1.0 || h <= hmin * 1.01
            t += h
            y .= ynew
            push!(raw_t, t)
            push!(raw_y, copy(y))
            dy_end = zeros(Float64, n)
            rhs!(dy_end, y, t)
            push!(raw_dy, dy_end)
            factor = err_norm == 0.0 ? 2.5 : min(2.5, 0.9 * err_norm^(-0.2))
            h = clamp(h * factor, hmin, hmax)
        else
            factor = max(0.2, 0.9 * err_norm^(-0.25))
            h = max(hmin, h * factor)
        end
    end

    if sample_times === nothing
        mat = zeros(Float64, n, length(raw_t))
        for k in eachindex(raw_t)
            mat[:, k] .= raw_y[k]
        end
        return raw_t, mat
    end

    ts = collect(Float64, sample_times)
    m = length(ts)
    mat = zeros(Float64, n, m)
    idx = 1
    for j in 1:m
        tj = clamp(ts[j], t0, tf)
        while idx < length(raw_t) - 1 && raw_t[idx + 1] < tj
            idx += 1
        end
        t_lo = raw_t[idx]
        t_hi = raw_t[min(idx + 1, length(raw_t))]
        dt_seg = t_hi - t_lo
        if dt_seg <= 1e-15
            mat[:, j] .= raw_y[idx]
        else
            # 4th-order cubic Hermite dense output interpolation
            s = (tj - t_lo) / dt_seg
            s2 = s * s
            s3 = s2 * s
            h00 = 2.0 * s3 - 3.0 * s2 + 1.0
            h10 = s3 - 2.0 * s2 + s
            h01 = -2.0 * s3 + 3.0 * s2
            h11 = s3 - s2
            @inbounds for i in 1:n
                mat[i, j] = h00 * raw_y[idx][i] +
                            h10 * dt_seg * raw_dy[idx][i] +
                            h01 * raw_y[idx + 1][i] +
                            h11 * dt_seg * raw_dy[idx + 1][i]
            end
        end
    end
    return ts, mat
end

"""
    finite_diff_jacobian!(J, rhs!, y, t; eps=1e-7)

Compute the finite-difference Jacobian `J = df/dy` at `(y, t)` in-place.
"""
function finite_diff_jacobian!(
    J::Matrix{Float64},
    rhs!::Function,
    y::Vector{Float64},
    t::Float64;
    eps::Float64 = 1e-7,
)
    n = length(y)
    f0 = zeros(Float64, n)
    f1 = zeros(Float64, n)
    yp = copy(y)
    rhs!(f0, y, t)
    for j in 1:n
        h = eps * max(1.0, abs(y[j]))
        yp[j] = y[j] + h
        rhs!(f1, yp, t)
        @inbounds for i in 1:n
            J[i, j] = (f1[i] - f0[i]) / h
        end
        yp[j] = y[j]
    end
    return J
end

"""
    implicit_euler_step(rhs!, y, t, dt; max_newton=15, tol=1e-9)

Advance `dy/dt = rhs!(dy, y, t)` by one stiff Backward Euler step `dt` using
damped Newton-Raphson with finite-difference Jacobian evaluation.
"""
function implicit_euler_step(
    rhs!::Function,
    y::Vector{Float64},
    t::Float64,
    dt::Float64;
    max_newton::Int = 15,
    tol::Float64 = 1e-9,
)
    n = length(y)
    t_next = t + dt
    yk = copy(y)
    fk = zeros(Float64, n)
    G = zeros(Float64, n)
    J = zeros(Float64, n, n)
    I_n = Matrix{Float64}(I, n, n)

    for _ in 1:max_newton
        rhs!(fk, yk, t_next)
        @inbounds for i in 1:n
            G[i] = yk[i] - y[i] - dt * fk[i]
        end
        if norm(G, Inf) <= tol
            return yk
        end
        finite_diff_jacobian!(J, rhs!, yk, t_next)
        A = I_n .- dt .* J
        delta = -(A \ G)
        yk .+= delta
        if norm(delta, Inf) <= tol
            return yk
        end
    end
    return yk
end

"""
    implicit_euler_integrate(rhs!, y0, tspan, n_steps)

Integrate a stiff ODE system across `tspan` with `n_steps` uniform Backward Euler steps.
"""
function implicit_euler_integrate(
    rhs!::Function,
    y0::AbstractVector{Float64},
    tspan::Tuple{Float64,Float64},
    n_steps::Int,
)
    n_steps >= 1 || throw(ArgumentError("n_steps must be >= 1"))
    t0, tf = Float64(tspan[1]), Float64(tspan[2])
    dt = (tf - t0) / n_steps
    n = length(y0)
    times = collect(range(t0, tf; length = n_steps + 1))
    states = zeros(Float64, n, n_steps + 1)
    y = copy(Float64.(y0))
    states[:, 1] .= y
    for k in 1:n_steps
        y = implicit_euler_step(rhs!, y, times[k], dt)
        states[:, k + 1] .= y
    end
    return times, states
end

"""
    crank_nicolson_integrate(rhs!, y0, tspan, n_steps; max_newton=15, tol=1e-9)

Second-order A-stable Crank-Nicolson trapezoidal implicit ODE integrator for
moderately stiff systems.
"""
function crank_nicolson_integrate(
    rhs!::Function,
    y0::AbstractVector{Float64},
    tspan::Tuple{Float64,Float64},
    n_steps::Int;
    max_newton::Int = 15,
    tol::Float64 = 1e-9,
)
    n_steps >= 1 || throw(ArgumentError("n_steps must be >= 1"))
    t0, tf = Float64(tspan[1]), Float64(tspan[2])
    dt = (tf - t0) / n_steps
    n = length(y0)
    times = collect(range(t0, tf; length = n_steps + 1))
    states = zeros(Float64, n, n_steps + 1)

    y = copy(Float64.(y0))
    states[:, 1] .= y
    f_prev = zeros(Float64, n)
    f_curr = zeros(Float64, n)
    G = zeros(Float64, n)
    J = zeros(Float64, n, n)
    I_n = Matrix{Float64}(I, n, n)

    for k in 1:n_steps
        t_k = times[k]
        t_next = times[k + 1]
        rhs!(f_prev, y, t_k)
        yk = y .+ dt .* f_prev
        for _ in 1:max_newton
            rhs!(f_curr, yk, t_next)
            @inbounds for i in 1:n
                G[i] = yk[i] - y[i] - 0.5 * dt * (f_prev[i] + f_curr[i])
            end
            if norm(G, Inf) <= tol
                break
            end
            finite_diff_jacobian!(J, rhs!, yk, t_next)
            A = I_n .- 0.5 .* dt .* J
            delta = -(A \ G)
            yk .+= delta
            if norm(delta, Inf) <= tol
                break
            end
        end
        y = yk
        states[:, k + 1] .= y
    end
    return times, states
end

"""
    solve_poisson_sor!(u, rhs, dx, dy; omega=1.72, tol=1e-6, max_iter=5000, dirichlet_zero=true)

Solve the 2D elliptic Poisson equation `d^2 u / dx^2 + d^2 u / dy^2 = rhs` on
a structured Cartesian grid using Successive Over-Relaxation (SOR).
Boundary points remain pinned if `dirichlet_zero=true`.
Returns `(iterations, max_residual)`.
"""
function solve_poisson_sor!(
    u::Matrix{Float64},
    rhs::Matrix{Float64},
    dx::Float64,
    dy::Float64;
    omega::Float64 = 1.72,
    tol::Float64 = 1e-6,
    max_iter::Int = 5000,
    dirichlet_zero::Bool = true,
)
    nx, ny = size(u)
    idx2 = 1.0 / (dx * dx)
    idy2 = 1.0 / (dy * dy)
    denom = 2.0 * (idx2 + idy2)

    if dirichlet_zero
        u[1, :] .= 0.0
        u[end, :] .= 0.0
        u[:, 1] .= 0.0
        u[:, end] .= 0.0
    end

    max_diff = Inf
    iters = 0
    for it in 1:max_iter
        iters = it
        max_diff = 0.0
        @inbounds for j in 2:(ny - 1)
            for i in 2:(nx - 1)
                gs = ((u[i - 1, j] + u[i + 1, j]) * idx2 +
                      (u[i, j - 1] + u[i, j + 1]) * idy2 - rhs[i, j]) / denom
                diff = omega * (gs - u[i, j])
                u[i, j] += diff
                ad = abs(diff)
                if ad > max_diff
                    max_diff = ad
                end
            end
        end
        if max_diff <= tol
            break
        end
    end
    return iters, max_diff
end

"""
    brent_root(f, a, b; tol=1e-10, max_iter=100)

Find a root of scalar function `f(x) = 0` bracketed in `[a, b]` using the
Brent-Dekker inverse quadratic / secant / bisection hybrid algorithm.
"""
function brent_root(f::Function, a::Float64, b::Float64; tol::Float64 = 1e-10, max_iter::Int = 100)
    fa = Float64(f(a))
    fb = Float64(f(b))
    abs(fa) <= tol && return a
    abs(fb) <= tol && return b
    fa * fb <= 0.0 || throw(ArgumentError("Root must be bracketed: f($a)=$fa, f($b)=$fb"))

    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end
    c = a
    fc = fa
    d = b - a
    e = d
    mflag = true

    for _ in 1:max_iter
        if abs(fb) <= tol || abs(b - a) <= tol
            return b
        end
        if fa != fc && fb != fc
            # Inverse quadratic interpolation
            s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
                (b * fa * fc) / ((fb - fa) * (fb - fc)) +
                (c * fa * fb) / ((fc - fa) * (fc - fb))
        else
            # Secant step
            s = b - fb * (b - a) / (fb - fa)
        end

        cond1 = !(min((3.0 * a + b) / 4.0, b) <= s <= max((3.0 * a + b) / 4.0, b))
        cond2 = mflag && abs(s - b) >= 0.5 * abs(b - c)
        cond3 = !mflag && abs(s - b) >= 0.5 * abs(c - d)
        cond4 = mflag && abs(b - c) < tol
        cond5 = !mflag && abs(c - d) < tol

        if cond1 || cond2 || cond3 || cond4 || cond5
            s = 0.5 * (a + b)
            mflag = true
        else
            mflag = false
        end

        fs = Float64(f(s))
        d = c
        c = b
        fc = fb
        if fa * fs < 0.0
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end
        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
    end
    return b
end

"""
    levenberg_marquardt(residual_fn, p0; lower=nothing, upper=nothing, max_iter=60, tol=1e-8)

Solve the nonlinear least-squares problem `min_p ||residual_fn(p)||_2^2` with
optional box bounds `[lower, upper]` using damped Levenberg-Marquardt updates.
Returns `(p_opt, rmse, iterations, converged)`.
"""
function levenberg_marquardt(
    residual_fn::Function,
    p0::AbstractVector{Float64};
    lower::Union{Nothing,AbstractVector{Float64}} = nothing,
    upper::Union{Nothing,AbstractVector{Float64}} = nothing,
    max_iter::Int = 60,
    tol::Float64 = 1e-8,
)
    p = copy(Float64.(p0))
    n = length(p)
    clamp_box!(v) = begin
        if lower !== nothing
            @inbounds for i in 1:n
                v[i] = max(v[i], lower[i])
            end
        end
        if upper !== nothing
            @inbounds for i in 1:n
                v[i] = min(v[i], upper[i])
            end
        end
        v
    end
    clamp_box!(p)

    r = Float64.(residual_fn(p))
    m = length(r)
    cost = 0.5 * dot(r, r)
    lambda = 1e-3
    converged = false
    iters = 0
    J = zeros(Float64, m, n)

    for it in 1:max_iter
        iters = it
        # Finite-difference Jacobian
        for j in 1:n
            h = 1e-6 * max(1.0, abs(p[j]))
            pj = copy(p)
            pj[j] += h
            clamp_box!(pj)
            dh = pj[j] - p[j]
            if abs(dh) < 1e-12
                pj = copy(p)
                pj[j] -= h
                clamp_box!(pj)
                dh = pj[j] - p[j]
            end
            rj = Float64.(residual_fn(pj))
            J[:, j] .= (rj .- r) ./ dh
        end

        JtJ = J' * J
        g = J' * r
        if norm(g, Inf) <= tol
            converged = true
            break
        end

        diag_scale = Diagonal(max.(diag(JtJ), 1e-6))
        accepted = false
        for _ in 1:10
            A = JtJ + lambda * diag_scale
            step = -(A \ g)
            p_cand = clamp_box!(p .+ step)
            r_cand = Float64.(residual_fn(p_cand))
            cost_cand = 0.5 * dot(r_cand, r_cand)
            if cost_cand < cost
                rel_impr = (cost - cost_cand) / max(1.0, cost)
                p = p_cand
                r = r_cand
                cost = cost_cand
                lambda = max(lambda * 0.3, 1e-9)
                accepted = true
                if rel_impr <= tol && norm(step, Inf) <= tol
                    converged = true
                end
                break
            else
                lambda = min(lambda * 5.0, 1e9)
            end
        end
        if converged || (!accepted && lambda >= 1e8)
            break
        end
    end
    rmse = sqrt(mean(r .^ 2))
    return p, rmse, iters, converged
end

"""
    WelfordAccumulator

Numerically stable single-pass streaming mean and sample variance accumulator.
"""
mutable struct WelfordAccumulator
    n::Int
    mean::Float64
    m2::Float64
    min_val::Float64
    max_val::Float64

    WelfordAccumulator() = new(0, 0.0, 0.0, Inf, -Inf)
end

function push_sample!(acc::WelfordAccumulator, x::Float64)
    acc.n += 1
    delta = x - acc.mean
    acc.mean += delta / acc.n
    delta2 = x - acc.mean
    acc.m2 += delta * delta2
    if x < acc.min_val
        acc.min_val = x
    end
    if x > acc.max_val
        acc.max_val = x
    end
    return acc
end

sample_variance(acc::WelfordAccumulator) = acc.n > 1 ? acc.m2 / (acc.n - 1) : 0.0
sample_std(acc::WelfordAccumulator) = sqrt(sample_variance(acc))

"""
    morris_sensitivity(model_fn, lower, upper; num_trajectories=12, levels=4, rng=MersenneTwister(42))

Perform Morris Elementary Effects global sensitivity analysis over box bounds
`[lower, upper]`. Returns `(mu_star, mu, sigma)` vectors quantifying mean absolute
influence and nonlinear/interaction spread for each input parameter.
"""
function morris_sensitivity(
    model_fn::Function,
    lower::Vector{Float64},
    upper::Vector{Float64};
    num_trajectories::Int = 12,
    levels::Int = 4,
    rng::AbstractRNG = MersenneTwister(42),
)
    d = length(lower)
    d == length(upper) || throw(ArgumentError("lower and upper bounds must match length"))
    delta_unit = levels / (2.0 * (levels - 1))
    effects = [Float64[] for _ in 1:d]

    for _ in 1:num_trajectories
        base = [rand(rng, 0:(levels - 1)) / (levels - 1) for _ in 1:d]
        order = randperm(rng, d)
        x_curr = lower .+ base .* (upper .- lower)
        y_curr = Float64(model_fn(x_curr))

        for j in order
            u_next = copy(base)
            if base[j] + delta_unit <= 1.0 + 1e-12
                u_next[j] = min(1.0, base[j] + delta_unit)
                step_u = u_next[j] - base[j]
            else
                u_next[j] = max(0.0, base[j] - delta_unit)
                step_u = u_next[j] - base[j]
            end
            x_next = lower .+ u_next .* (upper .- lower)
            y_next = Float64(model_fn(x_next))
            ee = (y_next - y_curr) / step_u
            push!(effects[j], ee)
            base = u_next
            x_curr = x_next
            y_curr = y_next
        end
    end

    mu_star = [mean(abs.(effects[j])) for j in 1:d]
    mu = [mean(effects[j]) for j in 1:d]
    sigma = [length(effects[j]) > 1 ? std(effects[j]) : 0.0 for j in 1:d]
    return mu_star, mu, sigma
end
