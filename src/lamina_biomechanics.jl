"""
Poroviscoelastic Lamina Cribrosa (LC) and Fung hyperelastic peripapillary
scleral shell biomechanics, translaminar pressure gradient (TLPG), and
axoplasmic retrograde transport blockade.

Models:
1. Thin-walled spherical Laplace hoop stress in the peripapillary sclera
   with Fung-type exponential stiffening:
   `sigma_hoop = P * R_globe / (2 * h_sclera)`
2. Scleral canal radial expansion `delta_R_canal` due to hoop strain, creating
   in-plane radial tensile stretch across the lamina cribrosa.
3. Axisymmetric Kirchhoff-Mindlin / circular-plate radial discretization of the
   lamina cribrosa under the net translaminar pressure load `delta_P = IOP - ICP`:
   `D_lc * nabla_r^4 w - N_r * nabla_r^2 w = delta_P`
   coupled with Kelvin-Voigt viscoelastic creep state `w_visc(t)`.
"""

const MMHG_TO_KPA = 0.1333223684

"""
    translaminar_pressure_gradient(iop_mmHg, icp_mmHg, lamina_thickness_mm)

Compute the Translaminar Pressure Gradient `TLPG = (IOP - ICP) / h_LC` in mmHg/mm.
Positive values push the lamina cribrosa posteriorly and impede retrograde BDNF
axoplasmic transport in retinal ganglion cell axons.
"""
function translaminar_pressure_gradient(iop_mmHg::Float64, icp_mmHg::Float64, lamina_thickness_mm::Float64)
    return (iop_mmHg - icp_mmHg) / max(lamina_thickness_mm, 0.05)
end

"""
    retrolaminar_icp(bio, hour_of_day)

Compute retrolaminar cerebrospinal fluid intracranial pressure (mmHg) accounting
for upright vs supine nocturnal postural shift.
"""
function retrolaminar_icp(bio::BiomechanicsParameters, hour_of_day::Float64)
    supine = is_supine_hour(hour_of_day)
    return bio.icp_upright_mmHg + supine * (bio.icp_supine_mmHg - bio.icp_upright_mmHg)
end

"""
    scleral_hoop_mechanics(anatomy, bio, iop_mmHg)

Compute peripapillary scleral hoop stress (kPa), tangent Fung modulus (kPa),
and scleral canal radial expansion (um) at a given `iop_mmHg`.
"""
function scleral_hoop_mechanics(anatomy::OcularAnatomy, bio::BiomechanicsParameters, iop_mmHg::Float64)
    p_kpa = max(iop_mmHg, 0.0) * MMHG_TO_KPA
    hoop_stress_kpa = p_kpa * anatomy.globe_radius_mm / (2.0 * anatomy.scleral_thickness_mm)
    # Fung nonlinear strain-stiffening inversion: sigma = (E0 / beta) * (exp(beta * eps) - 1)
    # => eps = (1 / beta) * log(1 + beta * sigma / E0)
    beta = max(bio.fung_stiffness_beta, 0.1)
    hoop_strain = (1.0 / beta) * log(1.0 + beta * hoop_stress_kpa / bio.scleral_modulus_kPa)
    canal_expansion_um = hoop_strain * anatomy.lamina_radius_mm * 1000.0
    return hoop_stress_kpa, hoop_strain, canal_expansion_um
end

"""
    solve_lamina_plate_profile(anatomy, bio, iop_mmHg, icp_mmHg; n_nodes=25, creep_factor=1.0)

Solve the axisymmetric radial boundary-value problem for the lamina cribrosa
posterior bowing displacement `w(r)` (um) across `n_nodes` radial stations
from `r = 0` (disc center) to `r = R_LC` (clamped/elastically supported at
the scleral canal rim), with geometric radial tension coupling from scleral
hoop stretch and Fung stiffness.

Returns `(r_mm, w_um, eps_r, shear_rz)`.
"""
function solve_lamina_plate_profile(
    anatomy::OcularAnatomy,
    bio::BiomechanicsParameters,
    iop_mmHg::Float64,
    icp_mmHg::Float64;
    n_nodes::Int = 25,
    creep_factor::Float64 = 1.0,
)
    n_nodes >= 9 || throw(ArgumentError("n_nodes must be >= 9"))
    R = anatomy.lamina_radius_mm
    h = anatomy.lamina_thickness_mm
    r = collect(range(0.0, R; length = n_nodes))
    dr = R / (n_nodes - 1)

    _, hoop_strain, _ = scleral_hoop_mechanics(anatomy, bio, iop_mmHg)
    delta_p_kpa = (iop_mmHg - icp_mmHg) * MMHG_TO_KPA

    # Effective flexural rigidity D = E * h^3 / (12 * (1 - nu^2)) [kPa * mm^3]
    nu = bio.lamina_poisson
    E_eff = bio.lamina_modulus_kPa / max(creep_factor, 0.25)
    D = E_eff * h^3 / (12.0 * (1.0 - nu^2))
    # Membrane in-plane radial tension T_r = E_eff * h * hoop_strain / (1 - nu) + collagen pre-tension
    G_lc = E_eff / (2.0 * (1.0 + nu))
    T_r = E_eff * h * hoop_strain / (1.0 - nu) + 0.035 * G_lc * h

    # Axisymmetric Reissner-Mindlin circular plate under uniform load q = delta_P:
    # Effective bending + transverse shear + membrane stiffness
    K_eff = T_r + 16.0 * D / (R^2) + 0.12 * G_lc * h

    A = zeros(Float64, n_nodes, n_nodes)
    b = fill(delta_p_kpa, n_nodes)

    # Symmetry at r = 0: limit r->0 of (w'' + w'/r) = 2 w'' => -4 K_eff (w_2 - w_1) / dr^2 = q
    A[1, 1] = 4.0 * K_eff / (dr^2)
    A[1, 2] = -4.0 * K_eff / (dr^2)

    for i in 2:(n_nodes - 1)
        ri = r[i]
        c_minus = K_eff * (1.0 / (dr^2) - 1.0 / (2.0 * ri * dr))
        c_diag = -2.0 * K_eff / (dr^2)
        c_plus = K_eff * (1.0 / (dr^2) + 1.0 / (2.0 * ri * dr))
        A[i, i - 1] = -c_minus
        A[i, i] = -c_diag
        A[i, i + 1] = -c_plus
    end

    # Clamped/pinned boundary at scleral insertion w(R) = 0
    A[n_nodes, n_nodes] = 1.0
    b[n_nodes] = 0.0

    w_mm = A \ b
    w_um = w_mm .* 1000.0

    # Compute local radial tensile strain (stretching + bending surface strain) and transverse shear
    eps_r = zeros(Float64, n_nodes)
    shear_rz = zeros(Float64, n_nodes)
    for i in 1:n_nodes
        dw_dr = if i == 1
            0.0
        elseif i == n_nodes
            (w_mm[n_nodes] - w_mm[n_nodes - 1]) / dr
        else
            (w_mm[i + 1] - w_mm[i - 1]) / (2.0 * dr)
        end
        d2w_dr2 = if i == 1
            2.0 * (w_mm[2] - w_mm[1]) / (dr^2)
        elseif i == n_nodes
            (w_mm[n_nodes] - 2.0 * w_mm[n_nodes - 1] + w_mm[n_nodes - 2]) / (dr^2)
        else
            (w_mm[i + 1] - 2.0 * w_mm[i] + w_mm[i - 1]) / (dr^2)
        end
        # Von Karman membrane + outer-fiber bending strain
        eps_r[i] = hoop_strain + 0.5 * dw_dr^2 + 0.5 * h * abs(d2w_dr2)
        # Shear strain peaks near scleral canal insertion r -> R
        shear_rz[i] = abs(dw_dr) * (0.4 + 0.6 * (r[i] / R))
    end
    return r, w_um, eps_r, shear_rz
end

"""
    evaluate_lamina_state(anatomy, bio, iop_mmHg, hour_of_day; creep_factor=1.0)

Compute the full `LaminaState` mechanical and axonal transport assessment at
instantaneous `iop_mmHg` and `hour_of_day`.
"""
function evaluate_lamina_state(
    anatomy::OcularAnatomy,
    bio::BiomechanicsParameters,
    iop_mmHg::Float64,
    hour_of_day::Float64;
    creep_factor::Float64 = 1.0,
)
    icp = retrolaminar_icp(bio, hour_of_day)
    tlpg = translaminar_pressure_gradient(iop_mmHg, icp, anatomy.lamina_thickness_mm)
    hoop_stress, _, canal_exp = scleral_hoop_mechanics(anatomy, bio, iop_mmHg)
    r, w_um, eps_r, shear_rz = solve_lamina_plate_profile(
        anatomy,
        bio,
        iop_mmHg,
        icp;
        creep_factor = creep_factor,
    )

    peak_bowing = w_um[1]
    # Effective von Mises-like axonal strain combining radial stretch and shear
    eff_strain = sqrt.(eps_r .^ 2 .+ 0.75 .* shear_rz .^ 2)
    mean_strain = mean(eff_strain)

    # Hill sigmoid for retrograde axoplasmic BDNF transport block
    ratio = mean_strain / max(bio.strain_injury_threshold, 1e-4)
    block_frac = (ratio^3) / (1.0 + ratio^3)

    return LaminaState(
        iop_mmHg,
        icp,
        tlpg,
        hoop_stress,
        canal_exp,
        peak_bowing,
        r,
        w_um,
        eps_r,
        shear_rz,
        mean_strain,
        block_frac,
    )
end
