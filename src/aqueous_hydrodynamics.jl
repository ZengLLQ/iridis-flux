"""
Aqueous humor hydrodynamics, Goldmann outflow balance, Schlemm's canal collapse,
postural episcleral venous pressure shifts, Friedenwald globe compliance, and
smart glaucoma drainage device (GDD) valve mechanics.

Governing ODE for intraocular pressure `P(t)` (mmHg):
    dP/dt = K * P * (F_in(t, P) - F_trab(t, P) - F_u(t) - F_valve(P))

where `1 / (K * P)` is the nonlinear Friedenwald ocular compliance `dV/dP`
(uL / mmHg), `F_in` is circadian ciliary active secretion minus pseudofacility,
`F_trab` is pressure-dependent trabecular outflow with Schlemm's canal collapse
and pupillary block modulation, `F_u` is uveoscleral outflow, and `F_valve` is
the optional drainage shunt flow.
"""

"""
    is_supine_hour(hour_of_day; sleep_start=23.0, sleep_end=7.0)

Return smooth posture indicator in `[0, 1]` representing supine sleep posture
with 30-minute physiological transition ramps at bedtime and awakening.
"""
function is_supine_hour(hour_of_day::Float64; sleep_start::Float64 = 23.0, sleep_end::Float64 = 7.0)
    h = mod(hour_of_day, 24.0)
    ramp = 0.5 # 30 min smooth sigmoid transition
    if h >= sleep_start
        return 0.5 * (1.0 + tanh((h - sleep_start) / (0.25 * ramp)))
    elseif h <= sleep_end
        return 0.5 * (1.0 - tanh((h - sleep_end) / (0.25 * ramp)))
    else
        # Transition down after wake or up before sleep
        return 0.0
    end
end

"""
    circadian_secretion_rate(aq, hour_of_day, iop; secretion_suppression=0.0)

Compute aqueous humor production rate `F_in` (uL/min) accounting for:
1. Diurnal circadian modulation (approx. 45% nocturnal suppression during sleep).
2. Ciliary pseudofacility (`-C_ps * (iop - 15.0)` ultrafiltration feedback).
3. Pharmacological suppression fraction in `[0, 0.9]`.
"""
function circadian_secretion_rate(
    aq::AqueousParameters,
    hour_of_day::Float64,
    iop::Float64;
    secretion_suppression::Float64 = 0.0,
)
    h = mod(hour_of_day, 24.0)
    # Cosine wave peaking near circadian_phase_shift_h (10:00) and troughing ~02:00
    phase = 2.0 * pi * (h - aq.circadian_phase_shift_h) / 24.0
    amp = 0.5 * aq.nocturnal_suppression
    circ_factor = (1.0 - amp) + amp * cos(phase)
    active_secretion = aq.base_secretion_uL_min * circ_factor * (1.0 - clamp(secretion_suppression, 0.0, 0.90))
    pseudo = aq.pseudofacility_uL_min_mmHg * max(iop - 15.0, -5.0)
    return max(0.05, active_secretion - pseudo)
end

"""
    trabecular_outflow_facility(aq, iop; facility_multiplier=1.0)

Compute effective trabecular meshwork outflow facility `C_trab` (uL/min/mmHg)
incorporating:
1. Schlemm's canal lumen collapse when `iop > P_crit` via exponential collapse:
   `exp(-alpha * max(iop - P_crit, 0))` (Moses & Brubaker model).
2. Pupillary block / angle closure occlusion `(1 - pupillary_block_factor)`.
3. Pharmacological relaxation multiplier (e.g., ROCK inhibitor netarsudil).
"""
function trabecular_outflow_facility(
    aq::AqueousParameters,
    iop::Float64;
    facility_multiplier::Float64 = 1.0,
)
    # Schlemm's canal collapse retains 30% collector-channel patency floor
    collapse = 0.30 + 0.70 * exp(-aq.schlemm_collapse_alpha * max(iop - aq.schlemm_collapse_pcrit_mmHg, 0.0))
    open_fraction = 1.0 - aq.pupillary_block_factor
    return aq.trabecular_facility_uL_min_mmHg * facility_multiplier * open_fraction * collapse
end

"""
    episcleral_venous_pressure(aq, hour_of_day; evp_reduction_mmHg=0.0)

Compute episcleral venous pressure `P_ev` (mmHg) including supine nocturnal
hydrostatic elevation and pharmacological venodilation reduction.
"""
function episcleral_venous_pressure(
    aq::AqueousParameters,
    hour_of_day::Float64;
    evp_reduction_mmHg::Float64 = 0.0,
)
    supine = is_supine_hour(hour_of_day)
    p_ev = aq.episcleral_venous_pressure_mmHg + supine * aq.supine_evp_delta_mmHg - evp_reduction_mmHg
    return max(4.0, p_ev)
end

"""
    trabecular_outflow_rate(aq, iop, hour_of_day; facility_multiplier=1.0, evp_reduction_mmHg=0.0)

Compute trabecular drainage flow `F_trab = C_trab(IOP) * (IOP - P_ev)` in uL/min.
Prevents non-physiological reflux below `P_ev` via diode valve check (`max(..., 0.0)`).
"""
function trabecular_outflow_rate(
    aq::AqueousParameters,
    iop::Float64,
    hour_of_day::Float64;
    facility_multiplier::Float64 = 1.0,
    evp_reduction_mmHg::Float64 = 0.0,
)
    c_eff = trabecular_outflow_facility(aq, iop; facility_multiplier = facility_multiplier)
    p_ev = episcleral_venous_pressure(aq, hour_of_day; evp_reduction_mmHg = evp_reduction_mmHg)
    return c_eff * max(iop - p_ev, 0.0)
end

"""
    uveoscleral_outflow_rate(aq; uveo_multiplier=1.0)

Compute pressure-insensitive uveoscleral drainage `F_u` (uL/min), modulated by
prostaglandin analog extracellular matrix remodeling.
"""
function uveoscleral_outflow_rate(aq::AqueousParameters; uveo_multiplier::Float64 = 1.0)
    return aq.uveoscleral_outflow_uL_min * max(0.1, uveo_multiplier)
end

"""
    friedenwald_compliance(aq, iop)

Evaluate instantaneous globe volume-pressure compliance `dV/dP = 1 / (K * P)`
in uL/mmHg using Friedenwald's exponential ocular rigidity law.
"""
function friedenwald_compliance(aq::AqueousParameters, iop::Float64)
    p_safe = max(iop, 3.0)
    return 1.0 / (aq.friedenwald_rigidity * p_safe)
end

"""
    valve_outflow_rate(valve, iop, conductance)

Compute flow (uL/min) through an implanted glaucoma drainage device (GDD) with
opening crack-pressure hysteresis and hypotony safety cutoff.
"""
function valve_outflow_rate(valve::DrainageValveConfig, iop::Float64, conductance::Float64)
    if !valve.enabled || iop <= valve.hypotony_guard_mmHg || iop <= valve.opening_pressure_mmHg
        return 0.0
    end
    # Smooth soft-opening sigmoid over 1 mmHg above opening pressure
    gate = 0.5 * (1.0 + tanh((iop - valve.opening_pressure_mmHg) / 0.6))
    return gate * conductance * (iop - valve.opening_pressure_mmHg)
end

"""
    aqueous_dpdt(aq, iop, hour_of_day; secretion_suppression=0.0, facility_multiplier=1.0, uveo_multiplier=1.0, evp_reduction_mmHg=0.0, valve=DrainageValveConfig(), valve_conductance=0.0)

Compute instantaneous rate of change `d(IOP)/dt` in **mmHg/hour** (converting
uL/min flows via factor `60.0 min/h`).
"""
function aqueous_dpdt(
    aq::AqueousParameters,
    iop::Float64,
    hour_of_day::Float64;
    secretion_suppression::Float64 = 0.0,
    facility_multiplier::Float64 = 1.0,
    uveo_multiplier::Float64 = 1.0,
    evp_reduction_mmHg::Float64 = 0.0,
    valve::DrainageValveConfig = DrainageValveConfig(),
    valve_conductance::Float64 = 0.0,
)
    f_in = circadian_secretion_rate(aq, hour_of_day, iop; secretion_suppression = secretion_suppression)
    f_trab = trabecular_outflow_rate(
        aq,
        iop,
        hour_of_day;
        facility_multiplier = facility_multiplier,
        evp_reduction_mmHg = evp_reduction_mmHg,
    )
    f_u = uveoscleral_outflow_rate(aq; uveo_multiplier = uveo_multiplier)
    f_v = valve_outflow_rate(valve, iop, valve_conductance)

    net_uL_per_h = 60.0 * (f_in - f_trab - f_u - f_v)
    compliance = friedenwald_compliance(aq, iop)
    return net_uL_per_h / compliance
end

"""
    steady_state_iop(aq; hour_of_day=10.0, secretion_suppression=0.0, facility_multiplier=1.0, uveo_multiplier=1.0, evp_reduction_mmHg=0.0)

Solve the nonlinear Goldmann equilibrium equation `F_in(P) - F_trab(P) - F_u = 0`
at fixed `hour_of_day` using Brent-Dekker root finding over `[4.0, 95.0]` mmHg.
"""
function steady_state_iop(
    aq::AqueousParameters;
    hour_of_day::Float64 = 10.0,
    secretion_suppression::Float64 = 0.0,
    facility_multiplier::Float64 = 1.0,
    uveo_multiplier::Float64 = 1.0,
    evp_reduction_mmHg::Float64 = 0.0,
)
    residual(p) = begin
        f_in = circadian_secretion_rate(aq, hour_of_day, p; secretion_suppression = secretion_suppression)
        f_trab = trabecular_outflow_rate(
            aq,
            p,
            hour_of_day;
            facility_multiplier = facility_multiplier,
            evp_reduction_mmHg = evp_reduction_mmHg,
        )
        f_u = uveoscleral_outflow_rate(aq; uveo_multiplier = uveo_multiplier)
        f_in - f_trab - f_u
    end
    # Scan for the lowest physiological sign-change bracket in [4.0, 140.0] mmHg
    p_grid = range(4.0, 140.0; length = 69)
    p_prev = first(p_grid)
    r_prev = residual(p_prev)
    for p_curr in p_grid[2:end]
        r_curr = residual(p_curr)
        if r_prev * r_curr <= 0.0
            return brent_root(residual, p_prev, p_curr; tol = 1e-8)
        end
        p_prev = p_curr
        r_prev = r_curr
    end
    return clamp(p_prev, 4.0, 140.0)
end

"""
    simulate_water_drinking_test(aq, baseline_iop; fluid_load_mL=1000.0, duration_min=60.0)

Simulate the clinical Water-Drinking Test (WDT) provocative stress challenge,
where rapid ingestion of 1 L hypotonic water induces transient episcleral venous
pressure elevation and ocular blood-volume expansion over 60 minutes, unmasking
impaired trabecular outflow facility. Returns `(time_min, iop_trace, peak_iop, iop_spike)`.
"""
function simulate_water_drinking_test(
    aq::AqueousParameters,
    baseline_iop::Float64;
    fluid_load_mL::Float64 = 1000.0,
    duration_min::Float64 = 60.0,
)
    scale = fluid_load_mL / 1000.0
    n_steps = Int(round(duration_min))
    dt_h = 1.0 / 60.0
    times_min = collect(0.0:1.0:duration_min)
    iops = zeros(Float64, length(times_min))
    iops[1] = baseline_iop

    p = baseline_iop
    for k in 1:n_steps
        t_m = times_min[k]
        # Transient osmotic choroidal expansion & EVP rise peaking ~20 min
        transient_evp_bump = 2.8 * scale * (t_m / 20.0) * exp(1.0 - t_m / 20.0)
        transient_inflow_bump = 0.35 * scale * exp(-t_m / 25.0)
        f_in = circadian_secretion_rate(aq, 10.0, p) + transient_inflow_bump
        c_eff = trabecular_outflow_facility(aq, p)
        p_ev = aq.episcleral_venous_pressure_mmHg + transient_evp_bump
        f_trab = c_eff * max(p - p_ev, 0.0)
        f_u = aq.uveoscleral_outflow_uL_min
        dpdt = 60.0 * (f_in - f_trab - f_u) / friedenwald_compliance(aq, p)
        p = clamp(p + dt_h * dpdt, 5.0, 80.0)
        iops[k + 1] = p
    end
    peak_iop = maximum(iops)
    return times_min, iops, peak_iop, peak_iop - baseline_iop
end
