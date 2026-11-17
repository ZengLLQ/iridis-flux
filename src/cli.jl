"""
Command-line interface (CLI) dispatcher for IridisFlux supporting 5 operational
subcommands: `simulate`, `anterior`, `calibrate`, `benchmark`, and `inspect`.
"""

using Printf

function _parse_cli_flags(args::Vector{String})
    opts = Dict{String,Any}(
        "phenotype" => "poag",
        "drug" => "none",
        "hours" => 24.0,
        "years" => 5.0,
        "steps" => 120,
        "valve" => false,
        "plot" => false,
        "help" => false,
        "target_iop" => 13.0,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            opts["help"] = true
        elseif arg == "--plot"
            opts["plot"] = true
        elseif arg == "--valve"
            opts["valve"] = true
        elseif arg == "--phenotype" && i < length(args)
            i += 1
            opts["phenotype"] = args[i]
        elseif arg == "--drug" && i < length(args)
            i += 1
            opts["drug"] = args[i]
        elseif arg == "--hours" && i < length(args)
            i += 1
            opts["hours"] = parse(Float64, args[i])
        elseif arg == "--years" && i < length(args)
            i += 1
            opts["years"] = parse(Float64, args[i])
        elseif arg == "--steps" && i < length(args)
            i += 1
            opts["steps"] = parse(Int, args[i])
        elseif arg == "--target-iop" && i < length(args)
            i += 1
            opts["target_iop"] = parse(Float64, args[i])
        end
        i += 1
    end
    return opts
end

function print_cli_help(io::IO = stdout)
    println(io, "IridisFlux - Multiscale Ocular Biophysics & Glaucoma Simulator")
    println(io, "")
    println(io, "Usage:")
    println(io, "  iridisflux <command> [options]")
    println(io, "")
    println(io, "Commands:")
    println(io, "  simulate   Run coupled diurnal IOP, PK/PD, lamina strain & 5-yr progression")
    println(io, "  anterior   Solve 2D sagittal Boussinesq anterior chamber convection PDE")
    println(io, "  calibrate  Fit patient trabecular facility & collapse from diurnal tonometry")
    println(io, "  benchmark  Run numerical solver convergence & performance benchmarks")
    println(io, "  inspect    Inspect clinical phenotype & therapy profiles (or full demo)")
    println(io, "")
    println(io, "Options:")
    println(io, "  --phenotype <name>   healthy | poag | ntg | aacc | pxfg (default: poag)")
    println(io, "  --drug <name>        none | latanoprost | timolol | dorzolamide | brimonidine | netarsudil | dual | max")
    println(io, "  --hours <float>      Diurnal simulation horizon in hours (default: 24.0)")
    println(io, "  --years <float>      Long-term progression horizon in years (default: 5.0)")
    println(io, "  --steps <int>        PDE time integration steps for anterior mode (default: 120)")
    println(io, "  --valve              Enable closed-loop PID MEMS drainage valve")
    println(io, "  --target-iop <float> Target IOP for closed-loop valve (default: 13.0)")
    println(io, "  --plot               Render 7-bit ASCII charts and perimetry maps")
    println(io, "  --help               Show this help message")
end

function run_cli(args::Vector{String}; io::IO = stdout)
    if isempty(args) || args[1] in ("-h", "--help")
        print_cli_help(io)
        return 0
    end

    cmd = lowercase(args[1])
    opts = _parse_cli_flags(args[2:end])
    if opts["help"]
        print_cli_help(io)
        return 0
    end

    if cmd == "simulate"
        patient = get_phenotype(opts["phenotype"])
        regimen = get_regimen(opts["drug"])
        valve = DrainageValveConfig(
            enabled = opts["valve"],
            closed_loop = opts["valve"],
            target_iop_mmHg = opts["target_iop"],
        )
        res = simulate_diurnal(patient; regimen = regimen, valve = valve, duration_h = opts["hours"])
        # Estimate aggregate neuroprotection from end-of-day aqueous PD
        neuro_bonus = isempty(regimen.agents) ? 0.0 : sum(a.neuroprotection_factor for a in regimen.agents) * 0.45
        prog = simulate_long_term_progression(patient, res; years = opts["years"], neuroprotection_bonus = neuro_bonus)
        risk = assess_glaucoma_risk(patient, res, prog)

        println(io, render_clinical_card(patient, res, prog, risk))
        if opts["plot"]
            println(io, "")
            println(io, render_iop_opp_chart(res))
            lam = evaluate_lamina_state(patient.anatomy, patient.biomechanics, res.peak_iop, 2.0)
            println(io, "")
            println(io, render_lamina_profile(lam))
            println(io, "")
            println(io, render_visual_field_24_2(prog.final_visual_field))
        end
        return 0

    elseif cmd == "anterior"
        patient = get_phenotype(opts["phenotype"])
        pde = solve_anterior_convection_2d(patient.anatomy; n_time_steps = opts["steps"], iop_mmHg = 24.0)
        @printf(io, "Anterior Chamber PDE Solved (%dx%d grid, %d steps)\n", pde.nx, pde.ny, pde.iterations)
        @printf(io, "  Peak Buoyancy Velocity  : %.4f mm/s\n", pde.max_velocity_mm_s)
        @printf(io, "  Corneal Nusselt Number  : %.3f\n", pde.mean_nusselt)
        @printf(io, "  Inferior Pigment Ratio  : %.2f%%\n", 100.0 * pde.inferior_deposition_ratio)
        @printf(io, "  Final Residual          : %.3e\n", pde.residual)
        if opts["plot"]
            println(io, "")
            println(io, render_anterior_chamber_map(pde))
        end
        return 0

    elseif cmd == "calibrate"
        # Synthetic nocturnal+diurnal tonometry profile from a hidden patient state
        true_patient = get_phenotype(opts["phenotype"])
        obs_hours = [8.0, 12.0, 16.0, 20.0, 2.0, 6.0]
        obs_iops = [steady_state_iop(true_patient.aqueous; hour_of_day = h) for h in obs_hours]
        _, fitted, rmse, iters, conv = calibrate_patient(obs_hours, obs_iops)
        println(io, "IridisFlux Levenberg-Marquardt Patient Outflow Calibration")
        @printf(io, "  Target Phenotype        : %s\n", true_patient.clinical_label)
        @printf(io, "  Converged               : %s (%d iterations, RMSE = %.4f mmHg)\n", string(conv), iters, rmse)
        @printf(io, "  Fitted Facility C0      : %.4f (True: %.4f uL/min/mmHg)\n", fitted.trabecular_facility, true_patient.aqueous.trabecular_facility_uL_min_mmHg)
        @printf(io, "  Fitted Uveoscleral Fu   : %.4f (True: %.4f uL/min)\n", fitted.uveoscleral_outflow, true_patient.aqueous.uveoscleral_outflow_uL_min)
        @printf(io, "  Fitted Collapse Alpha   : %.4f (True: %.4f 1/mmHg)\n", fitted.schlemm_collapse_alpha, true_patient.aqueous.schlemm_collapse_alpha)
        return 0

    elseif cmd == "benchmark"
        println(io, "IridisFlux Numerical Convergence & Performance Benchmark")
        # 1. RKF45 vs Analytical exponential decay order test
        f_exp!(dy, y, t) = (dy[1] = -0.35 * y[1]; dy)
        _, y_num = rkf45_integrate(f_exp!, [10.0], (0.0, 8.0); rtol = 1e-8, atol = 1e-10)
        err_rkf = abs(y_num[1, end] - 10.0 * exp(-0.35 * 8.0))

        # 2. Crank-Nicolson convergence test
        _, y_cn1 = crank_nicolson_integrate(f_exp!, [10.0], (0.0, 4.0), 40)
        _, y_cn2 = crank_nicolson_integrate(f_exp!, [10.0], (0.0, 4.0), 80)
        e1 = abs(y_cn1[1, end] - 10.0 * exp(-1.4))
        e2 = abs(y_cn2[1, end] - 10.0 * exp(-1.4))
        order_cn = log2(e1 / max(e2, 1e-16))

        # 3. Runtime benchmark on diurnal POAG simulation
        p = poag_moderate_phenotype()
        reg = get_regimen("latanoprost")
        t0 = time_ns()
        res = simulate_diurnal(p; regimen = reg, duration_h = 24.0)
        dt_ms = (time_ns() - t0) / 1e6

        # 4. Morris global sensitivity
        sens = run_aqueous_sensitivity(p; num_trajectories = 8)

        @printf(io, "  RKF45 Exponential Error : %.3e\n", err_rkf)
        @printf(io, "  Crank-Nicolson Order    : %.2f (expected ~2.00)\n", order_cn)
        @printf(io, "  24h Coupled Sim Time    : %.2f ms (Mean IOP = %.2f mmHg)\n", dt_ms, res.mean_iop)
        println(io, "  Morris Sensitivity Ranking (IOP drivers):")
        for row in sens
            @printf(io, "    %-25s mu* = %6.3f  sigma = %5.3f\n", row.parameter, row.mu_star, row.sigma)
        end
        return 0

    elseif cmd in ("inspect", "demo")
        patient = get_phenotype(opts["phenotype"])
        regimen = get_regimen(opts["drug"])
        res = simulate_diurnal(patient; regimen = regimen, duration_h = 24.0)
        prog = simulate_long_term_progression(patient, res; years = 5.0)
        risk = assess_glaucoma_risk(patient, res, prog)
        println(io, render_clinical_card(patient, res, prog, risk))
        if cmd == "demo" || opts["plot"]
            println(io, "")
            println(io, render_iop_opp_chart(res))
            println(io, "")
            println(io, render_visual_field_24_2(prog.final_visual_field))
        end
        return 0
    else
        println(io, "Unknown command: $cmd")
        print_cli_help(io)
        return 1
    end
end
