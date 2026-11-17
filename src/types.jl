"""
Core physical parameter structures, clinical phenotypes, and simulation outputs
for the IridisFlux multiscale ocular biophysics engine.

All physical quantities use consistent clinical and SI-derived units documented
explicitly on each field (pressures in mmHg, flows in uL/min, lengths in mm or um,
concentrations in ug/mL, time in hours or seconds as noted).
"""

"""
    OcularAnatomy

Macroscopic and microscopic geometric dimensions of the human eye globe,
anterior segment, and optic nerve head (ONH).

# Fields
- `globe_radius_mm::Float64`: Mean outer scleral shell radius (mm, normal ~12.0).
- `scleral_thickness_mm::Float64`: Peripapillary scleral thickness (mm, normal ~0.95).
- `anterior_chamber_depth_mm::Float64`: Central corneal endothelium to anterior lens pole (mm, normal ~3.1).
- `anterior_chamber_height_mm::Float64`: Vertical corneal limbus-to-limbus diameter (mm, normal ~11.5).
- `anterior_chamber_volume_uL::Float64`: Baseline anterior chamber fluid volume (uL, normal ~250.0).
- `corneal_thickness_mm::Float64`: Central corneal thickness (mm, normal ~0.545).
- `corneal_surface_area_cm2::Float64`: Exposed corneal epithelial absorption area (cm^2, normal ~1.04).
- `lamina_radius_mm::Float64`: Optic disc / scleral canal opening radius (mm, normal ~0.90).
- `lamina_thickness_mm::Float64`: Lamina cribrosa thickness across translaminar gradient (mm, normal ~0.25).
- `baseline_rnfl_um::Float64`: Healthy retinal nerve fiber layer thickness (um, normal ~102.0).
- `baseline_rgc_count::Float64`: Initial total retinal ganglion cell axon count (cells, normal ~1.1e6).
"""
struct OcularAnatomy
    globe_radius_mm::Float64
    scleral_thickness_mm::Float64
    anterior_chamber_depth_mm::Float64
    anterior_chamber_height_mm::Float64
    anterior_chamber_volume_uL::Float64
    corneal_thickness_mm::Float64
    corneal_surface_area_cm2::Float64
    lamina_radius_mm::Float64
    lamina_thickness_mm::Float64
    baseline_rnfl_um::Float64
    baseline_rgc_count::Float64

    function OcularAnatomy(;
        globe_radius_mm::Float64 = 12.0,
        scleral_thickness_mm::Float64 = 0.95,
        anterior_chamber_depth_mm::Float64 = 3.10,
        anterior_chamber_height_mm::Float64 = 11.5,
        anterior_chamber_volume_uL::Float64 = 250.0,
        corneal_thickness_mm::Float64 = 0.545,
        corneal_surface_area_cm2::Float64 = 1.04,
        lamina_radius_mm::Float64 = 0.90,
        lamina_thickness_mm::Float64 = 0.25,
        baseline_rnfl_um::Float64 = 102.0,
        baseline_rgc_count::Float64 = 1.10e6,
    )
        globe_radius_mm > 0.0 || throw(ArgumentError("globe_radius_mm must be positive"))
        scleral_thickness_mm > 0.0 || throw(ArgumentError("scleral_thickness_mm must be positive"))
        anterior_chamber_depth_mm > 0.0 || throw(ArgumentError("anterior_chamber_depth_mm must be positive"))
        anterior_chamber_volume_uL > 0.0 || throw(ArgumentError("anterior_chamber_volume_uL must be positive"))
        lamina_radius_mm > 0.0 || throw(ArgumentError("lamina_radius_mm must be positive"))
        lamina_thickness_mm > 0.0 || throw(ArgumentError("lamina_thickness_mm must be positive"))
        new(
            globe_radius_mm,
            scleral_thickness_mm,
            anterior_chamber_depth_mm,
            anterior_chamber_height_mm,
            anterior_chamber_volume_uL,
            corneal_thickness_mm,
            corneal_surface_area_cm2,
            lamina_radius_mm,
            lamina_thickness_mm,
            baseline_rnfl_um,
            baseline_rgc_count,
        )
    end
end

"""
    AqueousParameters

Parameters governing aqueous humor secretion, Schlemm's canal collapse,
uveoscleral drainage, episcleral venous pressure, and ocular globe rigidity.

# Fields
- `base_secretion_uL_min::Float64`: Daytime ciliary active secretion rate (uL/min, normal ~2.5).
- `nocturnal_suppression::Float64`: Fractional drop in secretion during sleep (0..0.8, normal ~0.45).
- `circadian_phase_shift_h::Float64`: Peak secretion phase offset in hours from midnight (normal ~10.0).
- `pseudofacility_uL_min_mmHg::Float64`: Pressure-dependent ciliary ultrafiltration suppression (uL/min/mmHg, normal ~0.02).
- `trabecular_facility_uL_min_mmHg::Float64`: Uncollapsed trabecular meshwork facility C0 (uL/min/mmHg, normal ~0.28).
- `schlemm_collapse_alpha::Float64`: Exponential Schlemm's canal collapse rate (1/mmHg, normal ~0.028).
- `schlemm_collapse_pcrit_mmHg::Float64`: Critical IOP threshold for canal narrowing (mmHg, normal ~18.0).
- `uveoscleral_outflow_uL_min::Float64`: Baseline pressure-insensitive uveoscleral drainage Fu (uL/min, normal ~0.45).
- `episcleral_venous_pressure_mmHg::Float64`: Upright daytime episcleral venous pressure Pev (mmHg, normal ~9.0).
- `supine_evp_delta_mmHg::Float64`: Postural hydrostatic EVP rise when supine during sleep (mmHg, normal ~3.2).
- `friedenwald_rigidity::Float64`: Friedenwald ocular rigidity coefficient K (1/uL, normal ~0.0215).
- `pupillary_block_factor::Float64`: Angle-closure trabecular occlusion multiplier in [0, 1) (0 = open angle, 0.92 = severe closure).
- `iris_lens_resistance::Float64`: Posterior-to-anterior chamber pupillary-block resistance (mmHg*min/uL).
"""
struct AqueousParameters
    base_secretion_uL_min::Float64
    nocturnal_suppression::Float64
    circadian_phase_shift_h::Float64
    pseudofacility_uL_min_mmHg::Float64
    trabecular_facility_uL_min_mmHg::Float64
    schlemm_collapse_alpha::Float64
    schlemm_collapse_pcrit_mmHg::Float64
    uveoscleral_outflow_uL_min::Float64
    episcleral_venous_pressure_mmHg::Float64
    supine_evp_delta_mmHg::Float64
    friedenwald_rigidity::Float64
    pupillary_block_factor::Float64
    iris_lens_resistance::Float64

    function AqueousParameters(;
        base_secretion_uL_min::Float64 = 2.50,
        nocturnal_suppression::Float64 = 0.45,
        circadian_phase_shift_h::Float64 = 10.0,
        pseudofacility_uL_min_mmHg::Float64 = 0.02,
        trabecular_facility_uL_min_mmHg::Float64 = 0.28,
        schlemm_collapse_alpha::Float64 = 0.028,
        schlemm_collapse_pcrit_mmHg::Float64 = 18.0,
        uveoscleral_outflow_uL_min::Float64 = 0.45,
        episcleral_venous_pressure_mmHg::Float64 = 9.0,
        supine_evp_delta_mmHg::Float64 = 3.2,
        friedenwald_rigidity::Float64 = 0.0215,
        pupillary_block_factor::Float64 = 0.0,
        iris_lens_resistance::Float64 = 0.35,
    )
        base_secretion_uL_min >= 0.0 || throw(ArgumentError("base_secretion_uL_min must be non-negative"))
        0.0 <= nocturnal_suppression <= 0.9 || throw(ArgumentError("nocturnal_suppression must be in [0, 0.9]"))
        trabecular_facility_uL_min_mmHg > 0.0 || throw(ArgumentError("trabecular_facility must be positive"))
        0.0 <= pupillary_block_factor < 1.0 || throw(ArgumentError("pupillary_block_factor must be in [0, 1)"))
        friedenwald_rigidity > 0.0 || throw(ArgumentError("friedenwald_rigidity must be positive"))
        new(
            base_secretion_uL_min,
            nocturnal_suppression,
            circadian_phase_shift_h,
            pseudofacility_uL_min_mmHg,
            trabecular_facility_uL_min_mmHg,
            schlemm_collapse_alpha,
            schlemm_collapse_pcrit_mmHg,
            uveoscleral_outflow_uL_min,
            episcleral_venous_pressure_mmHg,
            supine_evp_delta_mmHg,
            friedenwald_rigidity,
            pupillary_block_factor,
            iris_lens_resistance,
        )
    end
end

"""
    BiomechanicsParameters

Constitutive mechanical parameters for the peripapillary sclera and
poroviscoelastic lamina cribrosa (LC).

# Fields
- `scleral_modulus_kPa::Float64`: Tangent Young's modulus of peripapillary sclera (kPa, normal ~2800.0).
- `fung_stiffness_beta::Float64`: Fung exponential strain-stiffening exponent (dimensionless, normal ~4.5).
- `lamina_modulus_kPa::Float64`: Baseline lamina cribrosa shear/flexural modulus (kPa, normal ~220.0).
- `lamina_viscosity_kPa_h::Float64`: Kelvin-Voigt viscoelastic dashpot viscosity of LC (kPa*h, normal ~120.0).
- `lamina_poisson::Float64`: Effective drained Poisson ratio of LC collagen network (0..0.49, normal ~0.35).
- `icp_upright_mmHg::Float64`: Retrolaminar CSF intracranial pressure when upright (mmHg, normal ~7.5).
- `icp_supine_mmHg::Float64`: Retrolaminar CSF intracranial pressure when supine (mmHg, normal ~13.0).
- `remodeling_rate_per_year::Float64`: Collagen cross-linking / compliance remodeling rate (1/year).
- `strain_injury_threshold::Float64`: Critical axonal tensile/shear strain above which axoplasmic block accelerates.
"""
struct BiomechanicsParameters
    scleral_modulus_kPa::Float64
    fung_stiffness_beta::Float64
    lamina_modulus_kPa::Float64
    lamina_viscosity_kPa_h::Float64
    lamina_poisson::Float64
    icp_upright_mmHg::Float64
    icp_supine_mmHg::Float64
    remodeling_rate_per_year::Float64
    strain_injury_threshold::Float64

    function BiomechanicsParameters(;
        scleral_modulus_kPa::Float64 = 2800.0,
        fung_stiffness_beta::Float64 = 4.5,
        lamina_modulus_kPa::Float64 = 220.0,
        lamina_viscosity_kPa_h::Float64 = 120.0,
        lamina_poisson::Float64 = 0.35,
        icp_upright_mmHg::Float64 = 7.5,
        icp_supine_mmHg::Float64 = 13.0,
        remodeling_rate_per_year::Float64 = 0.04,
        strain_injury_threshold::Float64 = 0.045,
    )
        scleral_modulus_kPa > 0.0 || throw(ArgumentError("scleral_modulus_kPa must be positive"))
        lamina_modulus_kPa > 0.0 || throw(ArgumentError("lamina_modulus_kPa must be positive"))
        0.0 <= lamina_poisson < 0.5 || throw(ArgumentError("lamina_poisson must be in [0, 0.5)"))
        new(
            scleral_modulus_kPa,
            fung_stiffness_beta,
            lamina_modulus_kPa,
            lamina_viscosity_kPa_h,
            lamina_poisson,
            icp_upright_mmHg,
            icp_supine_mmHg,
            remodeling_rate_per_year,
            strain_injury_threshold,
        )
    end
end

"""
    RetinalVascularParameters

Parameters for retinal and optic nerve head (ONH) microvascular autoregulation
and Krogh tissue cylinder oxygen diffusion.

# Fields
- `map_day_mmHg::Float64`: Daytime brachial mean arterial pressure (mmHg, normal ~92.0).
- `map_night_mmHg::Float64`: Nocturnal brachial mean arterial pressure (mmHg, normal ~80.0; dipping < 70 increases NTG risk).
- `ophthalmic_map_ratio::Float64`: Ratio of ophthalmic artery pressure to brachial MAP upright (normal ~0.68).
- `supine_ophthalmic_ratio::Float64`: Ratio when supine without hydrostatic column loss (normal ~0.88).
- `autoregulation_gain::Float64`: Combined myogenic + metabolic + shear autoregulatory gain (0..2.5, normal ~1.0).
- `autoreg_tau_s::Float64`: Arteriolar vascular smooth muscle response time constant (seconds, normal ~18.0).
- `crit_opp_mmHg::Float64`: Lower limit of autoregulation ocular perfusion pressure (mmHg, normal ~32.0).
- `baseline_flow_uL_min::Float64`: Baseline ONH capillary perfusion flow (uL/min, normal ~18.5).
- `arterial_po2_mmHg::Float64`: Ophthalmic arteriolar plasma oxygen tension (mmHg, normal ~92.0).
- `krogh_radius_um::Float64`: Inter-capillary Krogh tissue cylinder outer radius (um, normal ~28.0).
- `capillary_radius_um::Float64`: ONH capillary lumen radius (um, normal ~3.2).
- `o2_consumption_max::Float64`: Michaelis-Menten Vmax tissue O2 consumption (mmHg/s equivalent, normal ~1.65).
"""
struct RetinalVascularParameters
    map_day_mmHg::Float64
    map_night_mmHg::Float64
    ophthalmic_map_ratio::Float64
    supine_ophthalmic_ratio::Float64
    autoregulation_gain::Float64
    autoreg_tau_s::Float64
    crit_opp_mmHg::Float64
    baseline_flow_uL_min::Float64
    arterial_po2_mmHg::Float64
    krogh_radius_um::Float64
    capillary_radius_um::Float64
    o2_consumption_max::Float64

    function RetinalVascularParameters(;
        map_day_mmHg::Float64 = 92.0,
        map_night_mmHg::Float64 = 80.0,
        ophthalmic_map_ratio::Float64 = 0.68,
        supine_ophthalmic_ratio::Float64 = 0.88,
        autoregulation_gain::Float64 = 1.0,
        autoreg_tau_s::Float64 = 18.0,
        crit_opp_mmHg::Float64 = 32.0,
        baseline_flow_uL_min::Float64 = 18.5,
        arterial_po2_mmHg::Float64 = 92.0,
        krogh_radius_um::Float64 = 28.0,
        capillary_radius_um::Float64 = 3.2,
        o2_consumption_max::Float64 = 1.65,
    )
        map_day_mmHg > 30.0 || throw(ArgumentError("map_day_mmHg must be > 30"))
        autoregulation_gain >= 0.0 || throw(ArgumentError("autoregulation_gain must be non-negative"))
        krogh_radius_um > capillary_radius_um || throw(ArgumentError("krogh_radius_um must exceed capillary_radius_um"))
        new(
            map_day_mmHg,
            map_night_mmHg,
            ophthalmic_map_ratio,
            supine_ophthalmic_ratio,
            autoregulation_gain,
            autoreg_tau_s,
            crit_opp_mmHg,
            baseline_flow_uL_min,
            arterial_po2_mmHg,
            krogh_radius_um,
            capillary_radius_um,
            o2_consumption_max,
        )
    end
end

"""
    DrugProfile

Mechanistic ocular pharmacokinetic and pharmacodynamic specification for a
topical or intracameral ophthalmic agent.
"""
struct DrugProfile
    name::String
    class_name::String
    dose_ug::Float64
    dosing_interval_h::Float64
    dosing_hours_of_day::Vector{Float64}
    corneal_permeability_per_h::Float64
    tear_drainage_per_h::Float64
    aqueous_clearance_per_h::Float64
    vitreous_transfer_per_h::Float64
    ec50_ug_mL::Float64
    hill_n::Float64
    emax_secretion_suppression::Float64
    emax_trabecular_boost::Float64
    emax_uveoscleral_boost::Float64
    emax_evp_reduction_mmHg::Float64
    neuroprotection_factor::Float64
    nocturnal_blunting::Float64
end

"""
    PharmacotherapyRegimen

Active multi-drug ocular regimen combining zero or more `DrugProfile` agents.
"""
struct PharmacotherapyRegimen
    name::String
    agents::Vector{DrugProfile}
    adherence_probability::Float64
end

"""
    DrainageValveConfig

Configuration for a passive or closed-loop MEMS glaucoma drainage device (GDD).
"""
struct DrainageValveConfig
    enabled::Bool
    closed_loop::Bool
    target_iop_mmHg::Float64
    opening_pressure_mmHg::Float64
    hypotony_guard_mmHg::Float64
    min_conductance_uL_min_mmHg::Float64
    max_conductance_uL_min_mmHg::Float64
    kp::Float64
    ki::Float64
    kd::Float64

    function DrainageValveConfig(;
        enabled::Bool = false,
        closed_loop::Bool = false,
        target_iop_mmHg::Float64 = 13.0,
        opening_pressure_mmHg::Float64 = 7.0,
        hypotony_guard_mmHg::Float64 = 6.0,
        min_conductance_uL_min_mmHg::Float64 = 0.02,
        max_conductance_uL_min_mmHg::Float64 = 0.65,
        kp::Float64 = 0.045,
        ki::Float64 = 0.015,
        kd::Float64 = 0.008,
    )
        new(
            enabled,
            closed_loop,
            target_iop_mmHg,
            opening_pressure_mmHg,
            hypotony_guard_mmHg,
            min_conductance_uL_min_mmHg,
            max_conductance_uL_min_mmHg,
            kp,
            ki,
            kd,
        )
    end
end

"""
    PatientPhenotype

Complete multi-scale patient parameter bundle characterizing ocular anatomy,
aqueous dynamics, lamina cribrosa biomechanics, and retinal perfusion.
"""
struct PatientPhenotype
    name::String
    clinical_label::String
    age_years::Float64
    anatomy::OcularAnatomy
    aqueous::AqueousParameters
    biomechanics::BiomechanicsParameters
    vascular::RetinalVascularParameters
end

"""
    AnteriorChamberPDEState

Output fields from the 2D sagittal Boussinesq anterior chamber thermal
convection and intracameral tracer transport PDE solver.
"""
struct AnteriorChamberPDEState
    nx::Int
    ny::Int
    width_mm::Float64
    height_mm::Float64
    x_mm::Vector{Float64}
    y_mm::Vector{Float64}
    temperature_C::Matrix{Float64}
    psi::Matrix{Float64}
    vorticity::Matrix{Float64}
    u_mm_s::Matrix{Float64}
    v_mm_s::Matrix{Float64}
    tracer_conc::Matrix{Float64}
    max_velocity_mm_s::Float64
    mean_nusselt::Float64
    inferior_deposition_ratio::Float64
    iterations::Int
    residual::Float64
end

"""
    LaminaState

Outputs of the axisymmetric poroviscoelastic lamina cribrosa & scleral shell
mechanical solver at a single time point.
"""
struct LaminaState
    iop_mmHg::Float64
    icp_mmHg::Float64
    tlpg_mmHg_mm::Float64
    scleral_hoop_stress_kPa::Float64
    scleral_canal_expansion_um::Float64
    peak_bowing_um::Float64
    radial_coords_mm::Vector{Float64}
    bowing_profile_um::Vector{Float64}
    radial_strain::Vector{Float64}
    shear_strain::Vector{Float64}
    mean_axonal_strain::Float64
    axoplasmic_block_fraction::Float64
end

"""
    PerfusionState

Outputs of the retinal/ONH microvascular autoregulation & Krogh O2 transport
model at a single time point.
"""
struct PerfusionState
    map_ophthalmic_mmHg::Float64
    opp_mmHg::Float64
    vascular_tone::Float64
    onh_blood_flow_uL_min::Float64
    flow_deficit_fraction::Float64
    mean_tissue_po2_mmHg::Float64
    min_tissue_po2_mmHg::Float64
    hypoxic_stress_index::Float64
end

"""
    VisualFieldMap

Humphrey 24-2 54-point perimetric sensitivity grid (dB) and global indices.
"""
struct VisualFieldMap
    grid_dB::Matrix{Float64}
    valid_mask::Matrix{Bool}
    sector_rgc_survival::Vector{Float64}
    sector_rnfl_um::Vector{Float64}
    mean_rnfl_um::Float64
    mean_deviation_dB::Float64
    pattern_std_dB::Float64
    visual_field_index_pct::Float64
end

"""
    DiurnalSimulationResult

Time-series output of coupled aqueous hydrodynamics, PK/PD, lamina biomechanics,
retinal perfusion, and smart valve control over a diurnal or multi-day window.
"""
struct DiurnalSimulationResult
    phenotype_name::String
    regimen_name::String
    time_h::Vector{Float64}
    iop_mmHg::Vector{Float64}
    icp_mmHg::Vector{Float64}
    opp_mmHg::Vector{Float64}
    secretion_uL_min::Vector{Float64}
    trabecular_outflow_uL_min::Vector{Float64}
    uveoscleral_outflow_uL_min::Vector{Float64}
    valve_outflow_uL_min::Vector{Float64}
    aqueous_drug_ug_mL::Vector{Float64}
    lamina_bowing_um::Vector{Float64}
    axonal_strain::Vector{Float64}
    onh_flow_uL_min::Vector{Float64}
    min_po2_mmHg::Vector{Float64}
    hypoxic_stress::Vector{Float64}
    valve_conductance::Vector{Float64}
    mean_iop::Float64
    peak_iop::Float64
    trough_iop::Float64
    iop_fluctuation::Float64
    min_opp::Float64
    peak_strain::Float64
    mean_hypoxic_stress::Float64
end

"""
    LongTermProgressionResult

Multi-year neurodegeneration forecast tracking RGC sector survival, RNFL
thinning, and Humphrey 24-2 perimetric indices.
"""
struct LongTermProgressionResult
    years::Vector{Float64}
    mean_iop_trajectory::Vector{Float64}
    mean_rnfl_um::Vector{Float64}
    mean_deviation_dB::Vector{Float64}
    vfi_pct::Vector{Float64}
    total_rgc_fraction::Vector{Float64}
    sector_survival_history::Matrix{Float64}
    final_visual_field::VisualFieldMap
    progression_rate_dB_yr::Float64
end
