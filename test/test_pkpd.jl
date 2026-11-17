@testset "Multi-Compartment Ocular PK/PD & Pharmacotherapy" begin
    anat = OcularAnatomy()
    poag = poag_moderate_phenotype()

    @testset "5-Compartment Tear-Cornea-Aqueous Mass Balance & Positivity" begin
        reg_lat = get_regimen("latanoprost")
        ts, aq_c, icb_c = simulate_ocular_pk(reg_lat, anat, 24.0; dt_sample_h = 0.25)
        @test length(ts) == 97
        @test minimum(aq_c) >= 0.0
        @test minimum(icb_c) >= 0.0
        @test maximum(aq_c) > 0.01
        @test maximum(icb_c) > 0.01

        # Verify dosing impulse integral equals single drop dose
        agent = latanoprost_agent()
        grid = collect(0.0:0.005:24.0)
        dose_integral = sum(dosing_impulse_rate(agent, h) for h in grid) * 0.005
        @test isapprox(dose_integral, agent.dose_ug; rtol = 0.03)
    end

    @testset "Efficacy Ordering & Nocturnal Timolol Blunting" begin
        res_none = simulate_diurnal(poag; regimen = no_therapy_regimen(), duration_h = 24.0)
        res_tim = simulate_diurnal(poag; regimen = get_regimen("timolol"), duration_h = 24.0)
        res_lat = simulate_diurnal(poag; regimen = get_regimen("latanoprost"), duration_h = 24.0)
        res_max = simulate_diurnal(poag; regimen = get_regimen("max"), duration_h = 24.0)

        @test res_tim.mean_iop < res_none.mean_iop - 2.0
        @test res_lat.mean_iop < res_none.mean_iop - 3.5
        @test res_max.mean_iop < res_lat.mean_iop

        # Timolol effect at 03:00 AM is blunted relative to 14:00 PM
        pd_day = combine_pharmacodynamics(get_regimen("timolol"), [3.0], 14.0)
        pd_night = combine_pharmacodynamics(get_regimen("timolol"), [3.0], 3.0)
        @test pd_day.secretion_suppression > 2.0 * pd_night.secretion_suppression
    end
end
