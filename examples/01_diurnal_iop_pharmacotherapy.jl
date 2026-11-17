"""
Example 01: Diurnal IOP Circadian Fluctuation & Glaucoma Pharmacotherapy Comparison

Compares an untreated Primary Open-Angle Glaucoma (POAG) patient across 48 hours
against Timolol 0.5% BID (nocturnal blunting), Latanoprost 0.005% QHS,
Netarsudil 0.02% QHS, and Maximal Triple Medical Therapy, plus a closed-loop
MEMS drainage valve.
"""

using Printf
using IridisFlux

function main()
    patient = poag_moderate_phenotype()
    regimen_names = ["none", "timolol", "latanoprost", "netarsudil", "max"]

    println("====================================================================")
    println("EXAMPLE 01: POAG DIURNAL IOP & MULTI-CLASS PHARMACOTHERAPY TRIAL")
    println("====================================================================")
    @printf("%-34s | %8s | %8s | %8s | %8s\n", "Regimen", "Mean IOP", "Peak IOP", "Min OPP", "5y MD")
    println(repeat("-", 76))

    last_res = nothing
    for rname in regimen_names
        reg = get_regimen(rname)
        res = simulate_diurnal(patient; regimen = reg, duration_h = 48.0, dt_sample_h = 0.25)
        neuro = isempty(reg.agents) ? 0.0 : sum(a.neuroprotection_factor for a in reg.agents) * 0.45
        prog = simulate_long_term_progression(patient, res; years = 5.0, neuroprotection_bonus = neuro)
        @printf(
            "%-34s | %6.2f   | %6.2f   | %6.2f   | %+6.2f dB\n",
            res.regimen_name[1:min(34, end)],
            res.mean_iop,
            res.peak_iop,
            res.min_opp,
            prog.mean_deviation_dB[end],
        )
        last_res = res
    end

    # Compare closed-loop MEMS drainage valve targeting 13.0 mmHg
    valve = DrainageValveConfig(enabled = true, closed_loop = true, target_iop_mmHg = 13.0)
    res_valve = simulate_diurnal(patient; valve = valve, duration_h = 48.0)
    prog_valve = simulate_long_term_progression(patient, res_valve; years = 5.0)
    @printf(
        "%-34s | %6.2f   | %6.2f   | %6.2f   | %+6.2f dB\n",
        "Closed-Loop MEMS Valve (13 mmHg)",
        res_valve.mean_iop,
        res_valve.peak_iop,
        res_valve.min_opp,
        prog_valve.mean_deviation_dB[end],
    )

    println("")
    if last_res !== nothing
        println(render_iop_opp_chart(last_res))
    end
end

main()
