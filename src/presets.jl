"""
Curated virtual clinical phenotypes and mechanistic ophthalmic drug profiles.
"""

"""
    healthy_emmetrope_phenotype()

Normal healthy adult eye (`IOP ~14-16 mmHg`, intact trabecular facility,
normal ICP, intact retinal microvascular autoregulation).
"""
function healthy_emmetrope_phenotype()
    return PatientPhenotype(
        "healthy",
        "Healthy Emmetrope (Normotensive Reference)",
        48.0,
        OcularAnatomy(),
        AqueousParameters(),
        BiomechanicsParameters(),
        RetinalVascularParameters(),
    )
end

"""
    poag_moderate_phenotype()

Primary Open-Angle Glaucoma (POAG) with markedly reduced trabecular outflow
facility (`C0 = 0.095 uL/min/mmHg`), increased Schlemm's canal collapsibility,
and compliant/aging lamina cribrosa (`IOP ~24-29 mmHg`).
"""
function poag_moderate_phenotype()
    return PatientPhenotype(
        "poag",
        "Primary Open-Angle Glaucoma (POAG)",
        67.0,
        OcularAnatomy(corneal_thickness_mm = 0.515, baseline_rnfl_um = 94.0),
        AqueousParameters(
            base_secretion_uL_min = 2.50,
            trabecular_facility_uL_min_mmHg = 0.135,
            schlemm_collapse_alpha = 0.025,
            schlemm_collapse_pcrit_mmHg = 16.0,
            uveoscleral_outflow_uL_min = 0.40,
            supine_evp_delta_mmHg = 3.8,
        ),
        BiomechanicsParameters(
            scleral_modulus_kPa = 2350.0,
            lamina_modulus_kPa = 175.0,
            icp_upright_mmHg = 7.0,
            strain_injury_threshold = 0.038,
        ),
        RetinalVascularParameters(
            map_day_mmHg = 90.0,
            map_night_mmHg = 76.0,
            autoregulation_gain = 0.72,
        ),
    )
end

"""
    normal_tension_glaucoma_phenotype()

Normal-Tension Glaucoma (NTG) with statistically normal IOP (~15-18 mmHg)
driven by low retrolaminar CSF intracranial pressure (`ICP = 3.2 mmHg` upright),
nocturnal systemic arterial hypotension (`MAP_night = 62 mmHg`), and impaired
autoregulation (`gain = 0.35`), causing high translaminar stress & ONH hypoxia.
"""
function normal_tension_glaucoma_phenotype()
    return PatientPhenotype(
        "ntg",
        "Normal-Tension Glaucoma (Low CSF + Nocturnal Dipper)",
        64.0,
        OcularAnatomy(lamina_thickness_mm = 0.19, baseline_rnfl_um = 92.0),
        AqueousParameters(
            base_secretion_uL_min = 2.35,
            trabecular_facility_uL_min_mmHg = 0.21,
            uveoscleral_outflow_uL_min = 0.42,
        ),
        BiomechanicsParameters(
            scleral_modulus_kPa = 2150.0,
            lamina_modulus_kPa = 155.0,
            icp_upright_mmHg = 3.2,
            icp_supine_mmHg = 7.5,
            strain_injury_threshold = 0.034,
        ),
        RetinalVascularParameters(
            map_day_mmHg = 84.0,
            map_night_mmHg = 62.0,
            autoregulation_gain = 0.35,
            crit_opp_mmHg = 35.0,
        ),
    )
end

"""
    acute_angle_closure_phenotype()

Acute Primary Angle-Closure Crisis (AACC) with shallow anterior chamber
(`2.05 mm`), severe trabecular meshwork occlusion (`pupillary_block_factor = 0.84`),
driving hypertensive spikes (`IOP > 45 mmHg`) and acute ischemic insult.
"""
function acute_angle_closure_phenotype()
    return PatientPhenotype(
        "aacc",
        "Acute Angle-Closure Crisis (Pupillary Block)",
        71.0,
        OcularAnatomy(anterior_chamber_depth_mm = 2.05, anterior_chamber_volume_uL = 165.0),
        AqueousParameters(
            base_secretion_uL_min = 2.60,
            trabecular_facility_uL_min_mmHg = 0.22,
            pupillary_block_factor = 0.84,
            uveoscleral_outflow_uL_min = 0.22,
        ),
        BiomechanicsParameters(
            scleral_modulus_kPa = 3100.0,
            lamina_modulus_kPa = 210.0,
        ),
        RetinalVascularParameters(
            map_day_mmHg = 96.0,
            map_night_mmHg = 84.0,
            autoregulation_gain = 0.65,
        ),
    )
end

"""
    pseudoexfoliation_phenotype()

Pseudoexfoliation Glaucoma (PXFG) characterized by flaky fibrillar debris
clogging the trabecular meshwork (`C0 = 0.075`), elevated diurnal IOP swing,
and weakened elastotic lamina cribrosa matrix.
"""
function pseudoexfoliation_phenotype()
    return PatientPhenotype(
        "pxfg",
        "Pseudoexfoliation Glaucoma (PXFG)",
        74.0,
        OcularAnatomy(baseline_rnfl_um = 90.0),
        AqueousParameters(
            base_secretion_uL_min = 2.55,
            nocturnal_suppression = 0.30,
            trabecular_facility_uL_min_mmHg = 0.105,
            schlemm_collapse_alpha = 0.032,
            uveoscleral_outflow_uL_min = 0.32,
            supine_evp_delta_mmHg = 4.1,
        ),
        BiomechanicsParameters(
            scleral_modulus_kPa = 1950.0,
            lamina_modulus_kPa = 140.0,
            strain_injury_threshold = 0.036,
        ),
        RetinalVascularParameters(
            map_day_mmHg = 90.0,
            map_night_mmHg = 75.0,
            autoregulation_gain = 0.55,
        ),
    )
end

function get_phenotype(name::AbstractString)
    key = lowercase(strip(name))
    if key in ("healthy", "normal", "emmetrope")
        return healthy_emmetrope_phenotype()
    elseif key in ("poag", "open_angle", "open-angle")
        return poag_moderate_phenotype()
    elseif key in ("ntg", "normal_tension", "normal-tension")
        return normal_tension_glaucoma_phenotype()
    elseif key in ("aacc", "angle_closure", "angle-closure")
        return acute_angle_closure_phenotype()
    elseif key in ("pxfg", "pseudoexfoliation")
        return pseudoexfoliation_phenotype()
    else
        throw(ArgumentError("Unknown phenotype '$name'. Valid options: healthy, poag, ntg, aacc, pxfg"))
    end
end

# ---------------------------------------------------------------------------
# Standard Ophthalmic Drugs & Combination Regimens
# ---------------------------------------------------------------------------

function latanoprost_agent()
    return DrugProfile(
        "Latanoprost 0.005%",
        "Prostaglandin F2alpha Analog",
        1.5,        # 30 uL of 50 ug/mL
        24.0,
        [21.0],     # QHS evening dosing
        0.55,
        8.0,
        0.52,
        0.04,
        0.035,
        1.2,
        0.0,        # no secretion suppression
        0.15,       # mild trabecular remodeling
        2.35,       # strong uveoscleral outflow enhancement
        0.0,
        0.05,
        0.0,        # active throughout night
    )
end

function timolol_agent()
    return DrugProfile(
        "Timolol 0.5%",
        "Non-Selective Beta-Blocker",
        150.0,      # 30 uL of 5000 ug/mL
        12.0,
        [8.0, 20.0],
        0.42,
        9.0,
        0.65,
        0.03,
        1.8,
        1.1,
        0.48,       # suppresses daytime ciliary secretion
        0.0,
        0.0,
        0.0,
        0.0,
        0.78,       # markedly blunted during nocturnal sleep
    )
end

function dorzolamide_agent()
    return DrugProfile(
        "Dorzolamide 2%",
        "Carbonic Anhydrase Inhibitor",
        600.0,
        8.0,
        [7.0, 15.0, 23.0],
        0.38,
        9.5,
        0.72,
        0.04,
        5.5,
        1.1,
        0.34,
        0.0,
        0.0,
        0.0,
        0.08,       # mild ONH vasodilatory CO2 effect
        0.10,
    )
end

function brimonidine_agent()
    return DrugProfile(
        "Brimonidine 0.2%",
        "Alpha-2 Adrenergic Agonist",
        60.0,
        12.0,
        [8.0, 20.0],
        0.48,
        8.5,
        0.68,
        0.05,
        0.75,
        1.15,
        0.30,
        0.0,
        0.45,
        0.0,
        0.42,       # direct BDNF / anti-apoptotic neuroprotection
        0.15,
    )
end

function netarsudil_agent()
    return DrugProfile(
        "Netarsudil 0.02%",
        "Rho-Kinase (ROCK) Inhibitor",
        6.0,
        24.0,
        [21.0],
        0.52,
        7.5,
        0.48,
        0.03,
        0.09,
        1.25,
        0.08,
        1.25,       # relaxes trabecular meshwork cytoskeleton
        0.0,
        2.1,        # lowers episcleral venous pressure
        0.10,
        0.05,
    )
end

no_therapy_regimen() = PharmacotherapyRegimen("Untreated Baseline", DrugProfile[], 1.0)

function get_regimen(name::AbstractString)
    key = lowercase(strip(name))
    if key in ("none", "untreated", "baseline")
        return no_therapy_regimen()
    elseif key in ("latanoprost", "pga")
        return PharmacotherapyRegimen("Latanoprost 0.005% QHS", [latanoprost_agent()], 0.92)
    elseif key in ("timolol", "beta", "betablocker")
        return PharmacotherapyRegimen("Timolol 0.5% BID", [timolol_agent()], 0.90)
    elseif key in ("dorzolamide", "cai")
        return PharmacotherapyRegimen("Dorzolamide 2% TID", [dorzolamide_agent()], 0.88)
    elseif key in ("brimonidine", "alpha2")
        return PharmacotherapyRegimen("Brimonidine 0.2% BID", [brimonidine_agent()], 0.88)
    elseif key in ("netarsudil", "rock")
        return PharmacotherapyRegimen("Netarsudil 0.02% QHS", [netarsudil_agent()], 0.92)
    elseif key in ("dual", "latanoprost+timolol", "combo")
        return PharmacotherapyRegimen("Latanoprost + Timolol Dual Therapy", [latanoprost_agent(), timolol_agent()], 0.90)
    elseif key in ("max", "triple", "maximal")
        return PharmacotherapyRegimen(
            "Maximal Medical Therapy (Latanoprost + Netarsudil + Brimonidine)",
            [latanoprost_agent(), netarsudil_agent(), brimonidine_agent()],
            0.88,
        )
    else
        throw(ArgumentError("Unknown regimen '$name'. Valid options: none, latanoprost, timolol, dorzolamide, brimonidine, netarsudil, dual, max"))
    end
end
