@testset "Neurodegeneration, Perimetry, Calibration & Valve Control" begin
    @testset "Humphrey 24-2 Visual Field & 5-Year Progression" begin
        healthy = healthy_emmetrope_phenotype()
        poag = poag_moderate_phenotype()

        vf_norm = compute_visual_field_24_2(ones(6), healthy.anatomy, 50.0)
        @test count(vf_norm.valid_mask) == 52
        @test isapprox(vf_norm.mean_deviation_dB, 0.0; atol = 0.2)
        @test vf_norm.visual_field_index_pct > 99.0

        res_h = simulate_diurnal(healthy; duration_h = 24.0)
        prog_h = simulate_long_term_progression(healthy, res_h; years = 5.0)
        res_p = simulate_diurnal(poag; duration_h = 24.0)
        prog_p = simulate_long_term_progression(poag, res_p; years = 5.0)

        @test prog_p.mean_deviation_dB[end] < prog_h.mean_deviation_dB[end] - 1.5
        @test prog_p.mean_rnfl_um[end] < prog_h.mean_rnfl_um[end]
        @test prog_p.vfi_pct[end] < prog_h.vfi_pct[end]
        @test prog_p.progression_rate_dB_yr < -0.2
    end

    @testset "Inverse Patient Outflow Calibration (Levenberg-Marquardt)" begin
        truth = poag_moderate_phenotype()
        obs_hours = [8.0, 12.0, 16.0, 20.0, 2.0, 6.0]
        obs_iops = [steady_state_iop(truth.aqueous; hour_of_day = h) for h in obs_hours]

        # Perturbed starting guess (same alpha/secretion prior, perturbed C0 and Fu)
        guess = PatientPhenotype(
            "guess",
            "Initial Guess",
            truth.age_years,
            truth.anatomy,
            AqueousParameters(
                base_secretion_uL_min = truth.aqueous.base_secretion_uL_min,
                trabecular_facility_uL_min_mmHg = 0.24,
                schlemm_collapse_alpha = truth.aqueous.schlemm_collapse_alpha,
                schlemm_collapse_pcrit_mmHg = truth.aqueous.schlemm_collapse_pcrit_mmHg,
                uveoscleral_outflow_uL_min = 0.70,
                supine_evp_delta_mmHg = truth.aqueous.supine_evp_delta_mmHg,
            ),
            truth.biomechanics,
            truth.vascular,
        )
        _, fitted, rmse, _, conv = calibrate_patient(obs_hours, obs_iops; template = guess)
        @test conv || rmse < 0.05
        @test rmse < 0.05
        @test isapprox(fitted.trabecular_facility, truth.aqueous.trabecular_facility_uL_min_mmHg; rtol = 0.05)
    end

    @testset "Closed-Loop PID MEMS Valve & 7-bit ASCII Guarantee" begin
        poag = poag_moderate_phenotype()
        valve = DrainageValveConfig(enabled = true, closed_loop = true, target_iop_mmHg = 13.0, hypotony_guard_mmHg = 6.0)
        res_v = simulate_diurnal(poag; valve = valve, duration_h = 48.0)

        @test res_v.trough_iop >= 6.0 # Hypotony guard never breached
        @test abs(res_v.iop_mmHg[end] - 13.0) < 2.5
        @test res_v.mean_iop < 17.0

        prog_v = simulate_long_term_progression(poag, res_v; years = 5.0)
        risk_v = assess_glaucoma_risk(poag, res_v, prog_v)
        lam = evaluate_lamina_state(poag.anatomy, poag.biomechanics, res_v.peak_iop, 12.0)
        pde = solve_anterior_convection_2d(poag.anatomy; nx = 17, ny = 25, n_time_steps = 40)

        # Verify strictly 7-bit ASCII across every visualizer output
        for txt in (
            render_iop_opp_chart(res_v),
            render_anterior_chamber_map(pde),
            render_lamina_profile(lam),
            render_visual_field_24_2(prog_v.final_visual_field),
            render_clinical_card(poag, res_v, prog_v, risk_v),
        )
            @test all(c -> UInt32(c) <= 127, txt)
        end
    end
end
