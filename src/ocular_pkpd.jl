"""
Multi-compartment ocular pharmacokinetics (PK) and mechanistic pharmacodynamics (PD)
for topical glaucoma therapeutics.

5-Compartment PK Topology (masses in ug, volumes in mL):
    1. Tear Film (`V_tear = 0.007 mL`): pulsed drop instillation + blinking/nasolacrimal washout
    2. Corneal Barrier (`V_cornea = 0.070 mL`): epithelium/stroma/endothelium reservoir
    3. Anterior Chamber Aqueous (`V_aq = 0.250 mL`): bulk turnover & trabecular clearance
    4. Iris-Ciliary Body (`V_icb = 0.150 mL`): target receptor site (ciliary epithelium, TM, uveoscleral)
    5. Posterior Vitreous (`V_vit = 4.000 mL`): slow deep-tissue sink

Pharmacodynamic coupling via Hill E_max equations with Bliss independence
multi-agent combination across prostaglandin analogs, beta-blockers, carbonic
anhydrase inhibitors, alpha-2 agonists, and ROCK inhibitors.
"""

const V_TEAR_ML = 0.007
const V_CORNEA_ML = 0.070
const V_ICB_ML = 0.150
const V_VITREOUS_ML = 4.000

"""
    dosing_impulse_rate(agent, hour_of_day; pulse_width_h=0.06)

Smooth Gaussian drop-instillation mass rate (ug/h) centered on scheduled
daily dosing hours (`agent.dosing_hours_of_day`), normalized so the integral
across each instillation equals `agent.dose_ug`.
"""
function dosing_impulse_rate(agent::DrugProfile, hour_of_day::Float64; pulse_width_h::Float64 = 0.05)
    h = mod(hour_of_day, 24.0)
    rate = 0.0
    norm_const = 1.0 / (sqrt(2.0 * pi) * pulse_width_h)
    for hd in agent.dosing_hours_of_day
        # Periodic distance on 24h clock
        dh = mod(h - hd + 12.0, 24.0) - 12.0
        rate += agent.dose_ug * norm_const * exp(-0.5 * (dh / pulse_width_h)^2)
    end
    return rate
end

"""
    pk_agent_rhs!(dy, y, t_h, agent, aq_volume_mL; adherence=1.0)

Right-hand side of the 5-compartment linear ocular PK ODE system for one `DrugProfile`:
`y = [m_tear, m_cornea, m_aqueous, m_icb, m_vitreous]` in micrograms (ug).
"""
function pk_agent_rhs!(
    dy::AbstractVector{Float64},
    y::AbstractVector{Float64},
    t_h::Float64,
    agent::DrugProfile,
    aq_volume_mL::Float64;
    adherence::Float64 = 1.0,
)
    m_tear, m_cornea, m_aq, m_icb, m_vit = y[1], y[2], y[3], y[4], y[5]

    in_rate = adherence * dosing_impulse_rate(agent, t_h)
    k_tc = agent.corneal_permeability_per_h
    k_drain = agent.tear_drainage_per_h
    k_ca = 1.65 * k_tc + 0.35 # stromal-to-aqueous release
    k_ac = 0.15 * k_tc        # back-diffusion aqueous->cornea
    k_ai = 0.95               # aqueous -> iris-ciliary body partition
    k_ia = 0.45               # iris-ciliary body -> aqueous washout
    k_clear = agent.aqueous_clearance_per_h
    k_av = agent.vitreous_transfer_per_h
    k_v_elim = 0.08           # retina-choroid systemic clearance from vitreous

    dy[1] = in_rate - (k_drain + k_tc) * m_tear
    dy[2] = k_tc * m_tear + k_ac * m_aq - (k_ca + 0.05) * m_cornea
    dy[3] = k_ca * m_cornea + k_ia * m_icb - (k_ac + k_ai + k_clear + k_av) * m_aq
    dy[4] = k_ai * m_aq - (k_ia + 0.12) * m_icb
    dy[5] = k_av * m_aq - k_v_elim * m_vit
    return dy
end

"""
    hill_effect(conc_ug_mL, ec50, hill_n)

Evaluate normalized Hill saturation response in `[0, 1]`.
"""
function hill_effect(conc_ug_mL::Float64, ec50::Float64, hill_n::Float64)
    c = max(conc_ug_mL, 0.0)
    c == 0.0 && return 0.0
    cn = c^hill_n
    return cn / (ec50^hill_n + cn)
end

"""
    combine_pharmacodynamics(regimen, icb_concentrations_ug_mL, hour_of_day)

Compute aggregate physiological pharmacodynamic modifiers across all active
agents in `regimen` at `hour_of_day`:
- `secretion_suppression`: fractional reduction in ciliary secretion (Bliss combined)
- `facility_multiplier`: multiplicative trabecular facility enhancement
- `uveo_multiplier`: multiplicative uveoscleral drainage enhancement
- `evp_reduction_mmHg`: additive episcleral venous pressure reduction (mmHg)
- `neuroprotection_bonus`: direct RGC axonal anti-apoptotic protection factor in `[0, 0.65]`
"""
function combine_pharmacodynamics(
    regimen::PharmacotherapyRegimen,
    icb_concentrations_ug_mL::AbstractVector{Float64},
    hour_of_day::Float64,
)
    supine = is_supine_hour(hour_of_day)
    surviving_secretion = 1.0
    facility_mult = 1.0
    uveo_mult = 1.0
    evp_red = 0.0
    surviving_neuro_risk = 1.0

    for (idx, agent) in enumerate(regimen.agents)
        c_icb = idx <= length(icb_concentrations_ug_mL) ? icb_concentrations_ug_mL[idx] : 0.0
        e_norm = hill_effect(c_icb, agent.ec50_ug_mL, agent.hill_n)
        # Nocturnal blunting (e.g., beta-blockers lose efficacy during low nocturnal sympathetic tone)
        circ_mod = 1.0 - supine * agent.nocturnal_blunting

        sec_eff = clamp(agent.emax_secretion_suppression * e_norm * circ_mod, 0.0, 0.85)
        surviving_secretion *= (1.0 - sec_eff)

        facility_mult *= (1.0 + agent.emax_trabecular_boost * e_norm)
        uveo_mult *= (1.0 + agent.emax_uveoscleral_boost * e_norm)
        evp_red += agent.emax_evp_reduction_mmHg * e_norm
        surviving_neuro_risk *= (1.0 - clamp(agent.neuroprotection_factor * e_norm, 0.0, 0.6))
    end

    secretion_suppression = clamp(1.0 - surviving_secretion, 0.0, 0.90)
    neuroprotection_bonus = clamp(1.0 - surviving_neuro_risk, 0.0, 0.65)
    return (
        secretion_suppression = secretion_suppression,
        facility_multiplier = facility_mult,
        uveo_multiplier = uveo_mult,
        evp_reduction_mmHg = evp_red,
        neuroprotection_bonus = neuroprotection_bonus,
    )
end

"""
    simulate_ocular_pk(regimen, anatomy, duration_h; dt_sample_h=0.25)

Integrate the 5-compartment ocular PK system for every agent in `regimen`
over `[0, duration_h]` using adaptive `rkf45_integrate` (initialized over a
48h pre-run warm-up for steady periodic dosing).
Returns `(times_h, aq_conc_total_ug_mL, icb_conc_matrix)`.
"""
function simulate_ocular_pk(
    regimen::PharmacotherapyRegimen,
    anatomy::OcularAnatomy,
    duration_h::Float64;
    dt_sample_h::Float64 = 0.25,
)
    sample_times = collect(0.0:dt_sample_h:duration_h)
    n_t = length(sample_times)
    n_agents = length(regimen.agents)
    aq_vol_ml = anatomy.anterior_chamber_volume_uL * 1e-3

    aq_total = zeros(Float64, n_t)
    icb_mat = zeros(Float64, n_agents, n_t)

    if n_agents == 0
        return sample_times, aq_total, icb_mat
    end

    warmup_h = 48.0
    shifted_samples = sample_times .+ warmup_h

    for (a_idx, agent) in enumerate(regimen.agents)
        y0 = zeros(Float64, 5)
        rhs!(dy, y, t) = pk_agent_rhs!(dy, y, t, agent, aq_vol_ml; adherence = regimen.adherence_probability)
        _, states = rkf45_integrate(
            rhs!,
            y0,
            (0.0, warmup_h + duration_h);
            rtol = 1e-6,
            atol = 1e-9,
            hmax = 0.04,
            sample_times = shifted_samples,
        )
        @inbounds for k in 1:n_t
            aq_total[k] += max(0.0, states[3, k]) / aq_vol_ml
            icb_mat[a_idx, k] = max(0.0, states[4, k]) / V_ICB_ML
        end
    end

    return sample_times, aq_total, icb_mat
end
