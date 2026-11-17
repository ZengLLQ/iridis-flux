"""
Retinal and Optic Nerve Head (ONH) microvascular autoregulation and
1D Krogh tissue cylinder oxygen diffusion-consumption model.

Models:
1. Ocular Perfusion Pressure `OPP(t) = MAP_ophthalmic(t) - IOP(t)` accounting for
   diurnal arterial blood pressure dipping and upright-vs-supine hydrostatic
   column correction.
2. Dynamic vascular smooth muscle tone `alpha_tone` driven by:
   - Myogenic Bayliss response (constriction on pressure rise, dilation on drop)
   - Metabolic CO2/lactate vasodilatory drive when flow falls below demand
   - Endothelial shear-mediated nitric oxide relaxation
   Yields the classic Lassen autoregulation plateau above `crit_opp_mmHg` and
   pressure-passive ischemic flow collapse below `crit_opp_mmHg`.
3. Radial Krogh cylinder O2 diffusion ODE/BVP around an ONH capillary with
   Michaelis-Menten mitochondrial consumption, yielding minimum perivascular
   `PO2` and Hypoxic Stress Index (`HSI`).
"""

"""
    ophthalmic_arterial_pressure(vasc, hour_of_day)

Compute mean ophthalmic artery pressure (mmHg) across the diurnal cycle,
combining circadian brachial MAP dipping with upright hydrostatic loss vs
supine perfusion gain.
"""
function ophthalmic_arterial_pressure(vasc::RetinalVascularParameters, hour_of_day::Float64)
    supine = is_supine_hour(hour_of_day)
    map_brachial = vasc.map_day_mmHg + supine * (vasc.map_night_mmHg - vasc.map_day_mmHg)
    ratio = vasc.ophthalmic_map_ratio + supine * (vasc.supine_ophthalmic_ratio - vasc.ophthalmic_map_ratio)
    return map_brachial * ratio
end

"""
    ocular_perfusion_pressure(vasc, iop_mmHg, hour_of_day)

Compute Ocular Perfusion Pressure `OPP = MAP_ophthalmic - IOP` in mmHg.
"""
function ocular_perfusion_pressure(vasc::RetinalVascularParameters, iop_mmHg::Float64, hour_of_day::Float64)
    p_oph = ophthalmic_arterial_pressure(vasc, hour_of_day)
    return p_oph - iop_mmHg
end

"""
    autoregulated_vascular_conductance(vasc, opp_mmHg; tone_override=nothing)

Compute steady-state or instantaneous normalized arteriolar conductance `G_norm`
in `[0.25, 2.2]` via Poiseuille `r^4` dilation/constriction around baseline
`OPP_0 = 46 mmHg`.

Above `crit_opp_mmHg`, autoregulation adjusts vascular tone so `Q = G_norm * (OPP / OPP_0)`
stays nearly flat on the Lassen plateau; below `crit_opp_mmHg`, maximal
vasodilation saturates and blood flow collapses linearly with `OPP`.
"""
function autoregulated_vascular_conductance(
    vasc::RetinalVascularParameters,
    opp_mmHg::Float64;
    baseline_opp::Float64 = 46.0,
)
    opp_safe = max(opp_mmHg, 1.0)
    ideal_conductance = baseline_opp / opp_safe
    # Max dilation limited by structural arteriolar lumen ceiling (~1.95x conductance)
    max_dilation = 1.0 + 0.95 * vasc.autoregulation_gain
    min_constriction = max(0.35, 1.0 - 0.55 * vasc.autoregulation_gain)

    # Smooth saturation near lower autoregulatory limit crit_opp_mmHg
    if opp_safe <= vasc.crit_opp_mmHg
        # Vessel approaches maximal vasodilation at crit_opp, then collapses
        collapse_factor = (opp_safe / vasc.crit_opp_mmHg)^0.35
        target = clamp(baseline_opp / vasc.crit_opp_mmHg, min_constriction, max_dilation)
        return target * collapse_factor
    else
        # Partial blend between constant flow and passive flow based on autoregulation_gain
        g_raw = 1.0 + vasc.autoregulation_gain * (ideal_conductance - 1.0)
        return clamp(g_raw, min_constriction, max_dilation)
    end
end

"""
    solve_krogh_oxygen_profile(vasc, blood_flow_uL_min; n_radial=24)

Solve the 1D axisymmetric Krogh tissue cylinder oxygen diffusion equation:
    (D_O2 * alpha_O2 / r) * d/dr( r * dP/dr ) = M(P)
where `M(P) = M_max * P / (K_m + P)` is Michaelis-Menten cellular respiration,
with convective capillary oxygen extraction reducing capillary wall PO2 when
`blood_flow_uL_min` drops.

Returns `(r_um, po2_profile_mmHg, mean_po2, min_po2)`.
"""
function solve_krogh_oxygen_profile(
    vasc::RetinalVascularParameters,
    blood_flow_uL_min::Float64;
    n_radial::Int = 24,
)
    n_radial >= 10 || throw(ArgumentError("n_radial must be >= 10"))
    r_cap = vasc.capillary_radius_um
    r_max = vasc.krogh_radius_um
    r = collect(range(r_cap, r_max; length = n_radial))
    dr = (r_max - r_cap) / (n_radial - 1)

    # Effective capillary inlet plasma PO2 drops with precapillary gradient and Fick extraction
    flow_ratio = clamp(blood_flow_uL_min / max(vasc.baseline_flow_uL_min, 1e-3), 0.05, 2.5)
    p_inlet = 0.72 * vasc.arterial_po2_mmHg
    extraction = 18.0 / flow_ratio
    p_cap = clamp(p_inlet - extraction, 5.0, p_inlet)

    # Effective Krogh diffusivity D_eff in (mmHg * um^2) / s giving realistic ~8-15 mmHg radial drop
    D_eff = 140.0
    km = 2.5

    # Solve nonlinear radial BVP via Picard / damped Newton iteration:
    # P(r_cap) = p_cap, dP/dr(r_max) = 0 (zero-flux symmetry between adjacent capillaries)
    p = fill(max(p_cap * 0.7, 10.0), n_radial)
    p[1] = p_cap

    A = zeros(Float64, n_radial, n_radial)
    rhs = zeros(Float64, n_radial)

    for _ in 1:25
        fill!(A, 0.0)
        fill!(rhs, 0.0)

        A[1, 1] = 1.0
        rhs[1] = p_cap

        for i in 2:(n_radial - 1)
            ri = r[i]
            # 1/r d/dr (r dP/dr) = d2P/dr2 + (1/r) dP/dr
            c_m = D_eff * (1.0 / dr^2 - 1.0 / (2.0 * ri * dr))
            c_0 = -2.0 * D_eff / dr^2
            c_p = D_eff * (1.0 / dr^2 + 1.0 / (2.0 * ri * dr))

            # Linearize Michaelis-Menten consumption M(P) = (M_max / (km + P_prev)) * P
            m_lin = vasc.o2_consumption_max / (km + max(p[i], 0.1))
            A[i, i - 1] = c_m
            A[i, i] = c_0 - m_lin
            A[i, i + 1] = c_p
            rhs[i] = 0.0
        end

        # Zero-flux Neumann condition at outer Krogh cylinder radius r = r_max
        A[n_radial, n_radial - 1] = -1.0
        A[n_radial, n_radial] = 1.0
        rhs[n_radial] = 0.0

        p_new = clamp.(A \ rhs, 0.05, p_cap)
        if maximum(abs.(p_new .- p)) < 1e-6
            p .= p_new
            break
        end
        p .= 0.5 .* p .+ 0.5 .* p_new
    end

    # Area-weighted mean PO2 across concentric radial rings
    weights = r .* dr
    mean_po2 = sum(p .* weights) / sum(weights)
    min_po2 = p[end]
    return r, p, mean_po2, min_po2
end

"""
    evaluate_perfusion_state(vasc, iop_mmHg, hour_of_day; tone_state=nothing)

Compute full `PerfusionState` coupling ophthalmic perfusion pressure,
autoregulated blood flow, and Krogh tissue oxygenation.
"""
function evaluate_perfusion_state(
    vasc::RetinalVascularParameters,
    iop_mmHg::Float64,
    hour_of_day::Float64;
    tone_state::Union{Nothing,Float64} = nothing,
)
    p_oph = ophthalmic_arterial_pressure(vasc, hour_of_day)
    opp = p_oph - iop_mmHg
    baseline_opp = 46.0

    g_eq = autoregulated_vascular_conductance(vasc, opp; baseline_opp = baseline_opp)
    g_eff = tone_state === nothing ? g_eq : clamp(tone_state, 0.25, 2.2)

    q_flow = max(0.2, vasc.baseline_flow_uL_min * g_eff * max(opp, 1.0) / baseline_opp)
    deficit = clamp(1.0 - q_flow / vasc.baseline_flow_uL_min, 0.0, 1.0)

    _, _, mean_po2, min_po2 = solve_krogh_oxygen_profile(vasc, q_flow)
    # Hypoxic Stress Index rises sharply as distal Krogh PO2 falls below 22 mmHg
    po2_crit = 22.0
    hsi = min_po2 >= po2_crit ? 0.0 : clamp(((po2_crit - min_po2) / po2_crit)^1.5, 0.0, 1.0)

    return PerfusionState(p_oph, opp, g_eff, q_flow, deficit, mean_po2, min_po2, hsi)
end
