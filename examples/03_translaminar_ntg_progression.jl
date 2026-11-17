"""
Example 03: Normal-Tension Glaucoma (NTG) Translaminar Pressure Gradient &
5-Year Humphrey 24-2 Visual Field Progression with Neuroprotective Rescue.

Contrasts a healthy normotensive eye (IOP ~15 mmHg, ICP 7.5 mmHg) against a
Normal-Tension Glaucoma phenotype (IOP ~16 mmHg, low CSF ICP 3.2 mmHg,
nocturnal systemic arterial dipping MAP_night = 62 mmHg), showing how
brimonidine + latanoprost suppresses both translaminar strain and hypoxia.
"""

using Printf
using IridisFlux

function main()
    println("====================================================================")
    println("EXAMPLE 03: NORMAL-TENSION GLAUCOMA & TRANSLAMINAR PROGRESSION")
    println("====================================================================")

    healthy = healthy_emmetrope_phenotype()
    ntg = normal_tension_glaucoma_phenotype()

    res_healthy = simulate_diurnal(healthy; duration_h = 24.0)
    prog_healthy = simulate_long_term_progression(healthy, res_healthy; years = 5.0)

    res_ntg = simulate_diurnal(ntg; duration_h = 24.0)
    prog_ntg = simulate_long_term_progression(ntg, res_ntg; years = 5.0)
    risk_ntg = assess_glaucoma_risk(ntg, res_ntg, prog_ntg)

    reg_rescue = get_regimen("max")
    res_rescue = simulate_diurnal(ntg; regimen = reg_rescue, duration_h = 24.0)
    prog_rescue = simulate_long_term_progression(ntg, res_rescue; years = 5.0, neuroprotection_bonus = 0.25)

    @printf("Healthy Reference : Mean IOP=%.2f mmHg | Min OPP=%.2f | 5y MD=%+.2f dB | VFI=%.1f%%\n",
        res_healthy.mean_iop, res_healthy.min_opp, prog_healthy.mean_deviation_dB[end], prog_healthy.vfi_pct[end])
    @printf("Untreated NTG     : Mean IOP=%.2f mmHg | Min OPP=%.2f | 5y MD=%+.2f dB | VFI=%.1f%%\n",
        res_ntg.mean_iop, res_ntg.min_opp, prog_ntg.mean_deviation_dB[end], prog_ntg.vfi_pct[end])
    @printf("Max Therapy NTG   : Mean IOP=%.2f mmHg | Min OPP=%.2f | 5y MD=%+.2f dB | VFI=%.1f%%\n\n",
        res_rescue.mean_iop, res_rescue.min_opp, prog_rescue.mean_deviation_dB[end], prog_rescue.vfi_pct[end])

    println(render_clinical_card(ntg, res_ntg, prog_ntg, risk_ntg))
    println("")
    println(render_visual_field_24_2(prog_ntg.final_visual_field))
end

main()
