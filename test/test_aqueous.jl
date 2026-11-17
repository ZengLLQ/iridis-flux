@testset "Aqueous Humor Hydrodynamics & Goldmann Balance" begin
    aq_healthy = AqueousParameters()
    aq_poag = poag_moderate_phenotype().aqueous
    aq_aacc = acute_angle_closure_phenotype().aqueous

    @testset "Steady-State Goldmann Conservation & Phenotype Ordering" begin
        p_h = steady_state_iop(aq_healthy; hour_of_day = 10.0)
        p_p = steady_state_iop(aq_poag; hour_of_day = 10.0)
        p_a = steady_state_iop(aq_aacc; hour_of_day = 10.0)

        @test 12.0 <= p_h <= 18.0
        @test 21.0 <= p_p <= 34.0
        @test p_a > 38.0
        @test p_h < p_p < p_a

        # Verify exact Goldmann flow conservation at equilibrium: F_in == F_trab + F_u
        f_in = circadian_secretion_rate(aq_healthy, 10.0, p_h)
        f_tr = trabecular_outflow_rate(aq_healthy, p_h, 10.0)
        f_u = uveoscleral_outflow_rate(aq_healthy)
        @test isapprox(f_in, f_tr + f_u; atol = 1e-6)
    end

    @testset "Schlemm's Canal Collapse & Supine EVP Shift" begin
        c_low = trabecular_outflow_facility(aq_healthy, 14.0)
        c_high = trabecular_outflow_facility(aq_healthy, 36.0)
        @test c_high < c_low

        p_ev_day = episcleral_venous_pressure(aq_healthy, 14.0)
        p_ev_night = episcleral_venous_pressure(aq_healthy, 2.0)
        @test p_ev_night > p_ev_day + 2.0
    end

    @testset "Water-Drinking Test Provocation & Valve Hypotony Guard" begin
        p0_h = steady_state_iop(aq_healthy; hour_of_day = 10.0)
        p0_p = steady_state_iop(aq_poag; hour_of_day = 10.0)
        _, _, _, spike_h = simulate_water_drinking_test(aq_healthy, p0_h)
        _, _, _, spike_p = simulate_water_drinking_test(aq_poag, p0_p)
        @test spike_h > 0.5
        @test spike_p > spike_h # Impaired trabecular facility causes larger WDT spike

        valve = DrainageValveConfig(enabled = true, opening_pressure_mmHg = 7.0, hypotony_guard_mmHg = 6.0)
        @test valve_outflow_rate(valve, 5.5, 0.4) == 0.0
        @test valve_outflow_rate(valve, 22.0, 0.4) > 3.0
    end
end
