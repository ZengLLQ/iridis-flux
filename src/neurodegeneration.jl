"""
Retinal Ganglion Cell (RGC) sectoral axonal neurodegeneration, RNFL thickness
structural decline, and Humphrey 24-2 54-point visual field perimetry mapping.

Maps coupled diurnal mechanical lamina cribrosa strain (`axoplasmic_block_fraction`)
and ischemic hypoxic stress (`hypoxic_stress_index`) onto the 6 Garway-Heath
anatomical sectors of the optic nerve head:
  1. Superotemporal (ST)  - highest biomechanical & arcuate vulnerability
  2. Inferotemporal (IT)  - highest early glaucomatous notch vulnerability
  3. Temporal / Maculopapillary (T)
  4. Superonasal (SN)
  5. Inferonasal (IN)
  6. Nasal (N)
"""

const GARWAY_HEATH_SECTORS = (
    "Superotemporal",
    "Inferotemporal",
    "Temporal",
    "Superonasal",
    "Inferonasal",
    "Nasal",
)

# Relative sector susceptibility weights reflecting clinical arcuate Bjerrum bundle vulnerability
const SECTOR_VULNERABILITY = Float64[1.28, 1.38, 0.62, 0.92, 0.96, 0.72]

# Standard 8x9 Humphrey 24-2 right-eye test pattern mask (54 valid visual field points)
const HUMPHREY_24_2_MASK = Bool[
    false false true  true  true  true  false false false;
    false true  true  true  true  true  true  false false;
    true  true  true  true  true  true  true  true  false;
    true  true  true  true  true  true  true  false true ;
    true  true  true  true  true  true  true  false true ;
    true  true  true  true  true  true  true  true  false;
    false true  true  true  true  true  true  false false;
    false false true  true  true  true  false false false;
]

# Garway-Heath anatomically mapped sector index (1..6) for each 24-2 grid location
# Upper hemifield maps to inferior disc (IT=2, IN=5); lower hemifield maps to superior disc (ST=1, SN=4)
const HUMPHREY_SECTOR_MAP = Int[
    0 0 2 2 2 5 0 0 0;
    0 2 2 2 2 5 5 0 0;
    2 2 2 2 3 5 5 6 0;
    2 2 3 3 3 6 6 0 6;
    1 1 3 3 3 6 6 0 6;
    1 1 1 1 3 4 4 6 0;
    0 1 1 1 1 4 4 0 0;
    0 0 1 1 1 4 0 0 0;
]

# Age-matched normal sensitivity hill-of-vision baseline (dB) at eccentricity (r, c)
function _normal_24_2_sensitivity(row::Int, col::Int, age_years::Float64)
    # Foveal-adjacent center at (4.5, 4.5); blind spot sits around (4:5, 8)
    ecc = sqrt((row - 4.5)^2 + (col - 4.5)^2)
    age_loss = 0.06 * max(age_years - 30.0, 0.0)
    return clamp(34.2 - 1.15 * ecc - age_loss, 18.0, 35.0)
end

"""
    compute_visual_field_24_2(sector_survival, anatomy, age_years)

Construct a `VisualFieldMap` (54-point Humphrey 24-2 dB grid, RNFL thickness,
Mean Deviation `MD`, Pattern Standard Deviation `PSD`, and Visual Field Index
`VFI %`) from 6-sector RGC survival fractions in `[0, 1]`.

Uses the Harwerth-Garway nonlinear log-linear structure-function relationship
where significant visual field loss (`dB`) accelerates once sector RGC survival
falls below ~75% ("pre-perimetric glaucoma" threshold).
"""
function compute_visual_field_24_2(
    sector_survival::AbstractVector{Float64},
    anatomy::OcularAnatomy,
    age_years::Float64,
)
    length(sector_survival) == 6 || throw(ArgumentError("sector_survival must have 6 entries"))
    surv = clamp.(Float64.(sector_survival), 0.01, 1.0)

    # OCT RNFL has a ~42 um non-axonal glial/vascular floor + axonal component
    glial_floor = 42.0
    axonal_span = max(anatomy.baseline_rnfl_um - glial_floor, 10.0)
    sector_rnfl = glial_floor .+ axonal_span .* surv
    mean_rnfl = mean(sector_rnfl)

    grid_dB = zeros(Float64, 8, 9)
    deviations = Float64[]

    for r in 1:8, c in 1:9
        if HUMPHREY_24_2_MASK[r, c]
            sec = HUMPHREY_SECTOR_MAP[r, c]
            s = surv[sec]
            norm_db = _normal_24_2_sensitivity(r, c, age_years)
            # Eccentricity-dependent local spatial gradient inside each arcuate zone
            ecc_mod = 1.0 + 0.05 * ((r - 4.5)^2 + (c - 4.5)^2) / 16.0
            # Pre-perimetric reserve up to ~20% loss, then logarithmic dB drop
            loss_frac = max(0.0, 1.0 - s)
            eff_loss = max(0.0, loss_frac - 0.06) / 0.94
            db_loss = clamp(32.0 * (eff_loss^1.35) * ecc_mod, 0.0, norm_db)
            val = max(0.0, norm_db - db_loss)
            grid_dB[r, c] = val
            push!(deviations, val - norm_db)
        end
    end

    md = mean(deviations)
    psd = std(deviations)
    # VFI: weighted percentage of surviving sensitivity normalized to 100%
    vfi = clamp(100.0 * (1.0 + md / 30.0), 0.0, 100.0)

    return VisualFieldMap(
        grid_dB,
        copy(HUMPHREY_24_2_MASK),
        surv,
        sector_rnfl,
        mean_rnfl,
        md,
        psd,
        vfi,
    )
end

"""
    simulate_long_term_progression(patient, diurnal_res; years=5.0, steps_per_year=12, neuroprotection_bonus=0.0)

Project multi-year RGC sectoral apoptosis, RNFL thinning, and Humphrey 24-2
visual field progression driven by the diurnal biomechanical (`peak_strain`,
`axoplasmic_block`) and ischemic (`mean_hypoxic_stress`) insult profile.
"""
function simulate_long_term_progression(
    patient::PatientPhenotype,
    diurnal_res::DiurnalSimulationResult;
    years::Float64 = 5.0,
    steps_per_year::Int = 12,
    neuroprotection_bonus::Float64 = 0.0,
)
    n_steps = max(1, Int(round(years * steps_per_year)))
    dt_yr = years / n_steps
    time_yr = collect(range(0.0, years; length = n_steps + 1))

    # Integrated daily biomechanical + ischemic insult rates
    mean_strain = mean(diurnal_res.axonal_strain)
    peak_strain = maximum(diurnal_res.axonal_strain)
    mean_hsi = diurnal_res.mean_hypoxic_stress

    # Baseline healthy physiological aging attrition ~0.4% per year
    base_aging_rate = 0.004
    # Pathological insult above threshold (healthy eyes stay below injury threshold)
    mech_excess = max(0.0, 0.5 * (mean_strain + peak_strain) - 0.85 * patient.biomechanics.strain_injury_threshold)
    mech_drive = 2.6 * (mech_excess / max(patient.biomechanics.strain_injury_threshold, 1e-4))^1.35
    ischemic_drive = 0.55 * mean_hsi

    prot_factor = 1.0 - clamp(neuroprotection_bonus, 0.0, 0.65)
    total_insult = prot_factor * (mech_drive + ischemic_drive)

    sector_surv = ones(Float64, 6)
    sector_hist = zeros(Float64, 6, n_steps + 1)
    sector_hist[:, 1] .= sector_surv

    iop_traj = fill(diurnal_res.mean_iop, n_steps + 1)
    rnfl_traj = zeros(Float64, n_steps + 1)
    md_traj = zeros(Float64, n_steps + 1)
    vfi_traj = zeros(Float64, n_steps + 1)
    rgc_traj = zeros(Float64, n_steps + 1)

    vf0 = compute_visual_field_24_2(sector_surv, patient.anatomy, patient.age_years)
    rnfl_traj[1] = vf0.mean_rnfl_um
    md_traj[1] = vf0.mean_deviation_dB
    vfi_traj[1] = vf0.visual_field_index_pct
    rgc_traj[1] = 1.0

    last_vf = vf0
    for k in 1:n_steps
        t_yr = time_yr[k + 1]
        # Progressive lamina cribrosa compliance stiffening/remodeling subtly alters strain
        remodel = 1.0 + patient.biomechanics.remodeling_rate_per_year * t_yr * 0.25
        for s in 1:6
            # Autocrine/paracrine TNF-alpha secondary degeneration feedback when survival < 0.85
            bystander = 1.0 + 0.45 * max(0.0, 0.85 - sector_surv[s])
            hazard = base_aging_rate + SECTOR_VULNERABILITY[s] * total_insult * remodel * bystander * 0.08
            sector_surv[s] = max(0.02, sector_surv[s] * exp(-hazard * dt_yr))
        end
        sector_hist[:, k + 1] .= sector_surv
        vf = compute_visual_field_24_2(sector_surv, patient.anatomy, patient.age_years + t_yr)
        rnfl_traj[k + 1] = vf.mean_rnfl_um
        md_traj[k + 1] = vf.mean_deviation_dB
        vfi_traj[k + 1] = vf.visual_field_index_pct
        rgc_traj[k + 1] = mean(sector_surv)
        last_vf = vf
    end

    md_rate = years > 0.0 ? (md_traj[end] - md_traj[1]) / years : 0.0

    return LongTermProgressionResult(
        time_yr,
        iop_traj,
        rnfl_traj,
        md_traj,
        vfi_traj,
        rgc_traj,
        sector_hist,
        last_vf,
        md_rate,
    )
end
