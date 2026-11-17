@testset "Retinal Microvascular Autoregulation & Krogh Oxygenation" begin
    vasc = RetinalVascularParameters()

    @testset "Lassen Autoregulation Plateau & Hypoperfusion Collapse" begin
        # Above crit_opp (32 mmHg), autoregulated flow varies weakly across IOP 12..22
        perf_14 = evaluate_perfusion_state(vasc, 14.0, 12.0)
        perf_22 = evaluate_perfusion_state(vasc, 22.0, 12.0)
        @test perf_14.opp_mmHg > 40.0
        @test abs(perf_14.onh_blood_flow_uL_min - perf_22.onh_blood_flow_uL_min) / perf_14.onh_blood_flow_uL_min < 0.08

        # Severe ocular hypertension pushes OPP below crit_opp -> ischemic drop
        perf_45 = evaluate_perfusion_state(vasc, 45.0, 12.0)
        @test perf_45.opp_mmHg < vasc.crit_opp_mmHg
        @test perf_45.onh_blood_flow_uL_min < 0.7 * perf_14.onh_blood_flow_uL_min
        @test perf_45.min_tissue_po2_mmHg < perf_14.min_tissue_po2_mmHg
        @test perf_45.hypoxic_stress_index > perf_14.hypoxic_stress_index
    end

    @testset "Krogh Radial Oxygen Diffusion Monotonicity" begin
        r, po2, mean_p, min_p = solve_krogh_oxygen_profile(vasc, 18.0)
        @test length(r) == length(po2)
        # PO2 decays monotonically from capillary wall (r_cap) to outer cylinder (r_max)
        @test all(diff(po2) .<= 1e-9)
        @test min_p == po2[end]
        @test min_p <= mean_p <= po2[1]

        _, _, _, min_p_low = solve_krogh_oxygen_profile(vasc, 5.5)
        @test min_p_low < min_p
    end
end
