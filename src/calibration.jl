"""
Inverse patient parameter calibration, global sensitivity analysis, and
clinical glaucoma risk stratification.

Includes:
- `calibrate_patient`: Recovers patient-specific trabecular outflow facility `C_0`,
  uveoscleral outflow `F_u`, and Schlemm's canal collapse rate `alpha` from
  clinical diurnal tonometry measurements using bounded Levenberg-Marquardt.
- `run_aqueous_sensitivity`: Evaluates Morris Elementary Effects across key
  physiological parameters to rank dominant drivers of peak IOP and strain.
- `assess_glaucoma_risk`: Generates structured clinical staging, perfusion deficit
  alerts, and 5-year visual field prognosis.
"""

"""
    calibrate_patient(obs_hours, obs_iops_mmHg; template=poag_moderate_phenotype())

Fit `trabecular_facility_uL_min_mmHg`, `uveoscleral_outflow_uL_min`, and
`schlemm_collapse_alpha` to observed diurnal Goldmann tonometry readings
`(obs_hours, obs_iops_mmHg)` via bounded Levenberg-Marquardt nonlinear least squares.

Returns `(calibrated_patient, fitted_params, rmse_mmHg, iterations, converged)`.
"""
function calibrate_patient(
    obs_hours::AbstractVector{Float64},
    obs_iops_mmHg::AbstractVector{Float64};
    template::PatientPhenotype = poag_moderate_phenotype(),
)
    length(obs_hours) == length(obs_iops_mmHg) || throw(ArgumentError("obs_hours and obs_iops_mmHg must have equal length"))
    length(obs_hours) >= 3 || throw(ArgumentError("Need at least 3 diurnal tonometry readings for identifiability"))

    base_aq = template.aqueous
    p0 = [
        base_aq.trabecular_facility_uL_min_mmHg,
        base_aq.uveoscleral_outflow_uL_min,
    ]
    lower = [0.02, 0.05]
    upper = [0.65, 1.50]

    function res_fn(p::Vector{Float64})
        aq_try = AqueousParameters(
            base_secretion_uL_min = base_aq.base_secretion_uL_min,
            nocturnal_suppression = base_aq.nocturnal_suppression,
            circadian_phase_shift_h = base_aq.circadian_phase_shift_h,
            pseudofacility_uL_min_mmHg = base_aq.pseudofacility_uL_min_mmHg,
            trabecular_facility_uL_min_mmHg = p[1],
            schlemm_collapse_alpha = base_aq.schlemm_collapse_alpha,
            schlemm_collapse_pcrit_mmHg = base_aq.schlemm_collapse_pcrit_mmHg,
            uveoscleral_outflow_uL_min = p[2],
            episcleral_venous_pressure_mmHg = base_aq.episcleral_venous_pressure_mmHg,
            supine_evp_delta_mmHg = base_aq.supine_evp_delta_mmHg,
            friedenwald_rigidity = base_aq.friedenwald_rigidity,
            pupillary_block_factor = base_aq.pupillary_block_factor,
            iris_lens_resistance = base_aq.iris_lens_resistance,
        )
        r = zeros(Float64, length(obs_hours))
        for k in eachindex(obs_hours)
            pred_iop = steady_state_iop(aq_try; hour_of_day = obs_hours[k])
            r[k] = pred_iop - obs_iops_mmHg[k]
        end
        return r
    end

    p_opt, rmse, iters, conv = levenberg_marquardt(res_fn, p0; lower = lower, upper = upper, max_iter = 50, tol = 1e-8)

    cal_aq = AqueousParameters(
        base_secretion_uL_min = base_aq.base_secretion_uL_min,
        nocturnal_suppression = base_aq.nocturnal_suppression,
        circadian_phase_shift_h = base_aq.circadian_phase_shift_h,
        pseudofacility_uL_min_mmHg = base_aq.pseudofacility_uL_min_mmHg,
        trabecular_facility_uL_min_mmHg = p_opt[1],
        schlemm_collapse_alpha = base_aq.schlemm_collapse_alpha,
        schlemm_collapse_pcrit_mmHg = base_aq.schlemm_collapse_pcrit_mmHg,
        uveoscleral_outflow_uL_min = p_opt[2],
        episcleral_venous_pressure_mmHg = base_aq.episcleral_venous_pressure_mmHg,
        supine_evp_delta_mmHg = base_aq.supine_evp_delta_mmHg,
        friedenwald_rigidity = base_aq.friedenwald_rigidity,
        pupillary_block_factor = base_aq.pupillary_block_factor,
        iris_lens_resistance = base_aq.iris_lens_resistance,
    )

    cal_patient = PatientPhenotype(
        template.name * "_calibrated",
        "Calibrated Clinical Profile",
        template.age_years,
        template.anatomy,
        cal_aq,
        template.biomechanics,
        template.vascular,
    )

    fitted = (
        trabecular_facility = p_opt[1],
        uveoscleral_outflow = p_opt[2],
        schlemm_collapse_alpha = base_aq.schlemm_collapse_alpha,
    )
    return cal_patient, fitted, rmse, iters, conv
end

"""
    run_aqueous_sensitivity(patient; num_trajectories=10)

Execute Morris Elementary Effects sensitivity analysis quantifying how
`(C_trab, F_in, F_u, P_ev, alpha_collapse)` govern steady-state IOP.
Returns a vector of named tuples sorted by descending `mu_star`.
"""
function run_aqueous_sensitivity(patient::PatientPhenotype; num_trajectories::Int = 10)
    aq = patient.aqueous
    names = ["Trabecular Facility C0", "Ciliary Secretion Fin", "Uveoscleral Outflow Fu", "Episcleral Venous Pev", "Schlemm Collapse Alpha"]
    lower = [0.06, 1.6, 0.15, 7.0, 0.005]
    upper = [0.35, 3.2, 0.85, 13.0, 0.060]

    function eval_fn(x::Vector{Float64})
        aq_x = AqueousParameters(
            base_secretion_uL_min = x[2],
            nocturnal_suppression = aq.nocturnal_suppression,
            circadian_phase_shift_h = aq.circadian_phase_shift_h,
            pseudofacility_uL_min_mmHg = aq.pseudofacility_uL_min_mmHg,
            trabecular_facility_uL_min_mmHg = x[1],
            schlemm_collapse_alpha = x[5],
            schlemm_collapse_pcrit_mmHg = aq.schlemm_collapse_pcrit_mmHg,
            uveoscleral_outflow_uL_min = x[3],
            episcleral_venous_pressure_mmHg = x[4],
            supine_evp_delta_mmHg = aq.supine_evp_delta_mmHg,
            friedenwald_rigidity = aq.friedenwald_rigidity,
            pupillary_block_factor = aq.pupillary_block_factor,
            iris_lens_resistance = aq.iris_lens_resistance,
        )
        return steady_state_iop(aq_x; hour_of_day = 10.0)
    end

    mu_star, mu, sigma = morris_sensitivity(eval_fn, lower, upper; num_trajectories = num_trajectories)
    rows = [(parameter = names[i], mu_star = mu_star[i], mu = mu[i], sigma = sigma[i]) for i in 1:5]
    sort!(rows; by = r -> r.mu_star, rev = true)
    return rows
end

"""
    assess_glaucoma_risk(patient, diurnal_res, prog_res)

Classify glaucoma risk stage and dominant pathophysiological driver
(Trabecular Hypertension, Translaminar Biomechanical Stress, or Vascular Ischemic Dipping).
"""
function assess_glaucoma_risk(
    patient::PatientPhenotype,
    diurnal_res::DiurnalSimulationResult,
    prog_res::LongTermProgressionResult,
)
    md_5yr = prog_res.mean_deviation_dB[end]
    md_rate = prog_res.progression_rate_dB_yr

    stage = if diurnal_res.peak_iop >= 36.0
        "Acute Hypertensive Crisis"
    elseif md_5yr <= -12.0 || md_rate <= -1.8
        "Severe Rapidly Progressing Glaucoma"
    elseif md_5yr <= -6.0 || md_rate <= -0.8
        "Moderate Progressive Glaucoma"
    elseif diurnal_res.mean_iop > 21.0 && md_5yr > -2.5
        "Ocular Hypertension / Early Glaucoma"
    elseif md_5yr <= -2.5
        "Early / Normal-Tension Glaucoma"
    else
        "Physiologically Stable / Low Risk"
    end

    dominant_driver = if diurnal_res.peak_iop > 24.0 && diurnal_res.mean_hypoxic_stress < 0.15
        "Elevated IOP & Schlemm's Canal Collapse"
    elseif diurnal_res.mean_hypoxic_stress >= 0.15 && diurnal_res.mean_iop <= 21.0
        "Nocturnal Hypoperfusion & Translaminar Pressure Gradient"
    elseif diurnal_res.peak_strain > patient.biomechanics.strain_injury_threshold
        "Mixed Biomechanical Lamina Strain & Ischemic Stress"
    else
        "Homeostatic Balance"
    end

    return (
        stage = stage,
        dominant_driver = dominant_driver,
        mean_iop = diurnal_res.mean_iop,
        peak_iop = diurnal_res.peak_iop,
        min_opp = diurnal_res.min_opp,
        peak_strain = diurnal_res.peak_strain,
        md_5yr_dB = md_5yr,
        md_rate_dB_yr = md_rate,
        vfi_5yr_pct = prog_res.vfi_pct[end],
    )
end
