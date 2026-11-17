"""
    IridisFlux

Multiscale computational ophthalmology and glaucoma biophysics engine coupling
circadian aqueous humor hydrodynamics, 2D sagittal Boussinesq anterior chamber
convection, poroviscoelastic lamina cribrosa biomechanics, retinal microvascular
autoregulation with Krogh tissue oxygenation, 5-compartment topical ocular PK/PD,
and Humphrey 24-2 retinal ganglion cell neurodegeneration forecasting.
"""
module IridisFlux

using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

include("types.jl")
include("numerics.jl")
include("aqueous_hydrodynamics.jl")
include("anterior_chamber_pde.jl")
include("lamina_biomechanics.jl")
include("retinal_perfusion.jl")
include("ocular_pkpd.jl")
include("neurodegeneration.jl")
include("presets.jl")
include("simulation.jl")
include("calibration.jl")
include("visualizer.jl")
include("cli.jl")

export OcularAnatomy,
       AqueousParameters,
       BiomechanicsParameters,
       RetinalVascularParameters,
       DrugProfile,
       PharmacotherapyRegimen,
       DrainageValveConfig,
       PatientPhenotype,
       AnteriorChamberPDEState,
       LaminaState,
       PerfusionState,
       VisualFieldMap,
       DiurnalSimulationResult,
       LongTermProgressionResult,
       # Numerics
       rkf45_integrate,
       implicit_euler_step,
       implicit_euler_integrate,
       crank_nicolson_integrate,
       solve_poisson_sor!,
       brent_root,
       levenberg_marquardt,
       WelfordAccumulator,
       push_sample!,
       sample_variance,
       sample_std,
       morris_sensitivity,
       # Aqueous hydrodynamics
       is_supine_hour,
       circadian_secretion_rate,
       trabecular_outflow_facility,
       episcleral_venous_pressure,
       trabecular_outflow_rate,
       uveoscleral_outflow_rate,
       friedenwald_compliance,
       valve_outflow_rate,
       aqueous_dpdt,
       steady_state_iop,
       simulate_water_drinking_test,
       # Anterior chamber PDE
       solve_anterior_convection_2d,
       # Lamina biomechanics
       translaminar_pressure_gradient,
       retrolaminar_icp,
       scleral_hoop_mechanics,
       solve_lamina_plate_profile,
       evaluate_lamina_state,
       # Retinal perfusion & Krogh O2
       ophthalmic_arterial_pressure,
       ocular_perfusion_pressure,
       autoregulated_vascular_conductance,
       solve_krogh_oxygen_profile,
       evaluate_perfusion_state,
       # Ocular PK/PD
       dosing_impulse_rate,
       pk_agent_rhs!,
       hill_effect,
       combine_pharmacodynamics,
       simulate_ocular_pk,
       # Neurodegeneration & Perimetry
       GARWAY_HEATH_SECTORS,
       compute_visual_field_24_2,
       simulate_long_term_progression,
       # Presets
       healthy_emmetrope_phenotype,
       poag_moderate_phenotype,
       normal_tension_glaucoma_phenotype,
       acute_angle_closure_phenotype,
       pseudoexfoliation_phenotype,
       get_phenotype,
       latanoprost_agent,
       timolol_agent,
       dorzolamide_agent,
       brimonidine_agent,
       netarsudil_agent,
       no_therapy_regimen,
       get_regimen,
       # Coupled simulation & Calibration
       simulate_diurnal,
       calibrate_patient,
       run_aqueous_sensitivity,
       assess_glaucoma_risk,
       # Visualizer & CLI
       render_iop_opp_chart,
       render_anterior_chamber_map,
       render_lamina_profile,
       render_visual_field_24_2,
       render_clinical_card,
       print_cli_help,
       run_cli

end # module
