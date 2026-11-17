@testset "2D Anterior Chamber Boussinesq Convection PDE" begin
    anat = OcularAnatomy()
    pde = solve_anterior_convection_2d(
        anat;
        nx = 25,
        ny = 37,
        t_cornea_C = 34.0,
        t_iris_C = 37.0,
        n_time_steps = 100,
        pigment_settling_mm_s = 0.03,
    )

    @test all(isfinite, pde.temperature_C)
    @test all(isfinite, pde.psi)
    @test all(isfinite, pde.tracer_conc)

    # Boundary temperatures enforced strictly
    @test all(isapprox.(pde.temperature_C[1, :], 34.0; atol = 1e-10))
    @test all(isapprox.(pde.temperature_C[end, :], 37.0; atol = 1e-10))

    # Thermal buoyancy creates nonzero convection plume and positive heat transport
    @test pde.max_velocity_mm_s > 0.01
    @test pde.mean_nusselt >= 0.9

    # Upward warm plume near iris (x -> W) and downward cool current near cornea (x -> 0)
    mid_j = div(pde.ny, 2)
    @test pde.v_mm_s[end - 2, mid_j] > 0.0
    @test pde.v_mm_s[3, mid_j] < 0.0

    # Pigment settling concentrates more tracer in inferior half (Krukenberg spindle)
    pde_nosettle = solve_anterior_convection_2d(
        anat;
        nx = 25,
        ny = 37,
        t_cornea_C = 34.0,
        t_iris_C = 37.0,
        n_time_steps = 100,
        pigment_settling_mm_s = 0.0,
    )
    @test pde.inferior_deposition_ratio > pde_nosettle.inferior_deposition_ratio
    @test 0.45 <= pde.inferior_deposition_ratio <= 1.0
    @test minimum(pde.tracer_conc) >= 0.0
end
