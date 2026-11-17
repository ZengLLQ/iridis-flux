"""
Coupled diurnal multiscale ocular simulator integrating:
1. Aqueous humor Goldmann-Schlemm ODE (`IOP(t)`)
2. Multi-compartment topical ocular PK/PD (`simulate_ocular_pk` + `combine_pharmacodynamics`)
3. Dynamic arteriolar autoregulatory vascular tone ODE (`d(tone)/dt`)
4. Kelvin-Voigt lamina cribrosa viscoelastic creep state (`d(creep)/dt`)
5. Closed-loop PID MEMS glaucoma drainage valve controller (`update_valve_conductance`)
"""

"""
    simulate_diurnal(patient; regimen=no_therapy_regimen(), valve=DrainageValveConfig(), duration_h=24.0, dt_sample_h=0.25, start_hour=8.0)

Run a full coupled multiscale diurnal simulation over `duration_h` hours starting
at `start_hour` (24h clock, default `8.0` = 08:00 AM).

Uses adaptive `rkf45_integrate` for the 5-compartment ocular PK system and
a tightly coupled 3-state biophysical ODE `[IOP, vascular_tone, lamina_creep]`
alongside discrete PID valve resistance updates with anti-windup and hypotony guard.
"""
function simulate_diurnal(
    patient::PatientPhenotype;
    regimen::PharmacotherapyRegimen = no_therapy_regimen(),
    valve::DrainageValveConfig = DrainageValveConfig(),
    duration_h::Float64 = 24.0,
    dt_sample_h::Float64 = 0.25,
    start_hour::Float64 = 8.0,
)
    duration_h > 0.0 || throw(ArgumentError("duration_h must be positive"))
    dt_sample_h > 0.0 || throw(ArgumentError("dt_sample_h must be positive"))

    # 1. Pre-integrate periodic steady-state 5-compartment ocular PK
    sample_rel_times = collect(0.0:dt_sample_h:duration_h)
    n_t = length(sample_rel_times)

    # Shift PK query onto clock hours [start_hour, start_hour + duration_h]
    # Simulate PK over [0, 48 + start_hour + duration_h] and sample at 48 + start_hour + rel_t
    n_agents = length(regimen.agents)
    aq_vol_ml = patient.anatomy.anterior_chamber_volume_uL * 1e-3
    aq_drug_conc = zeros(Float64, n_t)
    icb_conc_mat = zeros(Float64, n_agents, n_t)

    if n_agents > 0
        warmup = 48.0
        query_times = sample_rel_times .+ (warmup + start_hour)
        for (a_idx, agent) in enumerate(regimen.agents)
            y0_pk = zeros(Float64, 5)
            rhs_pk!(dy, y, t) = pk_agent_rhs!(dy, y, t, agent, aq_vol_ml; adherence = regimen.adherence_probability)
            _, states_pk = rkf45_integrate(
                rhs_pk!,
                y0_pk,
                (0.0, warmup + start_hour + duration_h);
                rtol = 1e-6,
                atol = 1e-9,
                hmax = 0.04,
                sample_times = query_times,
            )
            @inbounds for k in 1:n_t
                aq_drug_conc[k] += max(0.0, states_pk[3, k]) / aq_vol_ml
                icb_conc_mat[a_idx, k] = max(0.0, states_pk[4, k]) / V_ICB_ML
            end
        end
    end

    # Helper to interpolate ICB concentrations at arbitrary clock hour t_abs
    function pd_at(t_rel::Float64)
        if n_agents == 0
            return (
                secretion_suppression = 0.0,
                facility_multiplier = 1.0,
                uveo_multiplier = 1.0,
                evp_reduction_mmHg = 0.0,
                neuroprotection_bonus = 0.0,
            )
        end
        idx_f = clamp(t_rel / dt_sample_h + 1.0, 1.0, Float64(n_t))
        i0 = floor(Int, idx_f)
        i1 = min(i0 + 1, n_t)
        w = idx_f - i0
        c_vec = [(1.0 - w) * icb_conc_mat[a, i0] + w * icb_conc_mat[a, i1] for a in 1:n_agents]
        return combine_pharmacodynamics(regimen, c_vec, start_hour + t_rel)
    end

    # 2. Initialize IOP at untreated/treated starting equilibrium + 24h periodic settle
    pd0 = pd_at(0.0)
    iop_init = steady_state_iop(
        patient.aqueous;
        hour_of_day = start_hour,
        secretion_suppression = pd0.secretion_suppression,
        facility_multiplier = pd0.facility_multiplier,
        uveo_multiplier = pd0.uveo_multiplier,
        evp_reduction_mmHg = pd0.evp_reduction_mmHg,
    )
    opp_init = ocular_perfusion_pressure(patient.vascular, iop_init, start_hour)
    tone_init = autoregulated_vascular_conductance(patient.vascular, opp_init)
    creep_init = 1.0

    # 3. Integrate coupled biophysical system step-by-step with discrete PID valve updates
    iop_vec = zeros(Float64, n_t)
    icp_vec = zeros(Float64, n_t)
    opp_vec = zeros(Float64, n_t)
    sec_vec = zeros(Float64, n_t)
    trab_vec = zeros(Float64, n_t)
    uveo_vec = zeros(Float64, n_t)
    valve_flow_vec = zeros(Float64, n_t)
    bowing_vec = zeros(Float64, n_t)
    strain_vec = zeros(Float64, n_t)
    flow_vec = zeros(Float64, n_t)
    po2_vec = zeros(Float64, n_t)
    hsi_vec = zeros(Float64, n_t)
    cond_vec = zeros(Float64, n_t)

    # Settle over 24h periodic warmup if valve disabled, else start directly from iop_init
    state = [iop_init, tone_init, creep_init]
    pid_integral = 0.0
    prev_error = 0.0
    curr_cond = valve.enabled ? 0.5 * (valve.min_conductance_uL_min_mmHg + valve.max_conductance_uL_min_mmHg) : 0.0

    # Viscoelastic Kelvin-Voigt relaxation time constant (hours)
    tau_creep_h = max(patient.biomechanics.lamina_viscosity_kPa_h / patient.biomechanics.lamina_modulus_kPa, 0.25)
    tau_autoreg_h = max(patient.vascular.autoreg_tau_s / 3600.0, 0.002)

    function bio_rhs!(dy, y, t_rel, cond_val)
        p_iop = clamp(y[1], 3.0, 95.0)
        tone = clamp(y[2], 0.25, 2.2)
        creep = clamp(y[3], 0.5, 1.8)
        clock_h = start_hour + t_rel
        pd = pd_at(t_rel)

        dy[1] = aqueous_dpdt(
            patient.aqueous,
            p_iop,
            clock_h;
            secretion_suppression = pd.secretion_suppression,
            facility_multiplier = pd.facility_multiplier,
            uveo_multiplier = pd.uveo_multiplier,
            evp_reduction_mmHg = pd.evp_reduction_mmHg,
            valve = valve,
            valve_conductance = cond_val,
        )
        opp_now = ocular_perfusion_pressure(patient.vascular, p_iop, clock_h)
        tone_eq = autoregulated_vascular_conductance(patient.vascular, opp_now)
        dy[2] = (tone_eq - tone) / tau_autoreg_h

        # Lamina creep relaxes slowly toward modest strain-dependent compliance drift
        icp_now = retrolaminar_icp(patient.biomechanics, clock_h)
        tlpg_now = max(0.0, p_iop - icp_now)
        creep_target = 1.0 + 0.12 * clamp((tlpg_now - 8.0) / 25.0, 0.0, 1.0)
        dy[3] = (creep_target - creep) / tau_creep_h
        return dy
    end

    # Short 24h periodic settling pass so initial state lies on limit cycle
    for k in 1:Int(round(24.0 / dt_sample_h))
        t0 = (k - 1) * dt_sample_h
        t1 = k * dt_sample_h
        # Use implicit Backward Euler on short step for rock-solid stiff tone stability
        step_rhs!(dy, y, t) = bio_rhs!(dy, y, t, curr_cond)
        state = implicit_euler_step(step_rhs!, state, t0, t1 - t0)
    end

    # Main recorded simulation pass
    for k in 1:n_t
        t_rel = sample_rel_times[k]
        clock_h = start_hour + t_rel
        p_iop = clamp(state[1], 3.0, 95.0)
        tone = clamp(state[2], 0.25, 2.2)
        creep = clamp(state[3], 0.5, 1.8)
        pd = pd_at(t_rel)

        iop_vec[k] = p_iop
        cond_vec[k] = curr_cond
        sec_vec[k] = circadian_secretion_rate(
            patient.aqueous,
            clock_h,
            p_iop;
            secretion_suppression = pd.secretion_suppression,
        )
        trab_vec[k] = trabecular_outflow_rate(
            patient.aqueous,
            p_iop,
            clock_h;
            facility_multiplier = pd.facility_multiplier,
            evp_reduction_mmHg = pd.evp_reduction_mmHg,
        )
        uveo_vec[k] = uveoscleral_outflow_rate(patient.aqueous; uveo_multiplier = pd.uveo_multiplier)
        valve_flow_vec[k] = valve_outflow_rate(valve, p_iop, curr_cond)

        lam = evaluate_lamina_state(patient.anatomy, patient.biomechanics, p_iop, clock_h; creep_factor = creep)
        perf = evaluate_perfusion_state(patient.vascular, p_iop, clock_h; tone_state = tone)

        icp_vec[k] = lam.icp_mmHg
        bowing_vec[k] = lam.peak_bowing_um
        strain_vec[k] = lam.mean_axonal_strain
        opp_vec[k] = perf.opp_mmHg
        flow_vec[k] = perf.onh_blood_flow_uL_min
        po2_vec[k] = perf.min_tissue_po2_mmHg
        hsi_vec[k] = perf.hypoxic_stress_index

        if k < n_t
            # Closed-loop PID drainage valve update with anti-windup & hypotony cutoff
            if valve.enabled && valve.closed_loop
                if p_iop <= valve.hypotony_guard_mmHg + 0.5
                    curr_cond = valve.min_conductance_uL_min_mmHg
                    pid_integral *= 0.5
                else
                    err = p_iop - valve.target_iop_mmHg
                    d_err = (err - prev_error) / dt_sample_h
                    cand_int = clamp(pid_integral + err * dt_sample_h, -8.0, 8.0)
                    u_pid = 0.25 + valve.kp * err + valve.ki * cand_int + valve.kd * d_err
                    clamped = clamp(u_pid, valve.min_conductance_uL_min_mmHg, valve.max_conductance_uL_min_mmHg)
                    # Conditional integration anti-windup
                    if clamped == u_pid
                        pid_integral = cand_int
                    end
                    curr_cond = clamped
                    prev_error = err
                end
            elseif valve.enabled
                curr_cond = 0.5 * (valve.min_conductance_uL_min_mmHg + valve.max_conductance_uL_min_mmHg)
            end

            step_rhs2!(dy, y, t) = bio_rhs!(dy, y, t, curr_cond)
            state = implicit_euler_step(step_rhs2!, state, t_rel, dt_sample_h)
        end
    end

    mean_iop = mean(iop_vec)
    peak_iop = maximum(iop_vec)
    trough_iop = minimum(iop_vec)

    return DiurnalSimulationResult(
        patient.name,
        regimen.name,
        sample_rel_times,
        iop_vec,
        icp_vec,
        opp_vec,
        sec_vec,
        trab_vec,
        uveo_vec,
        valve_flow_vec,
        aq_drug_conc,
        bowing_vec,
        strain_vec,
        flow_vec,
        po2_vec,
        hsi_vec,
        cond_vec,
        mean_iop,
        peak_iop,
        trough_iop,
        peak_iop - trough_iop,
        minimum(opp_vec),
        maximum(strain_vec),
        mean(hsi_vec),
    )
end
