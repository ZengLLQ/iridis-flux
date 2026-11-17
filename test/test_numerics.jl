@testset "Numerical Solvers & Optimization" begin
    @testset "RKF45 Exponential & Harmonic Oscillator Convergence" begin
        rhs_exp!(dy, y, t) = (dy[1] = -0.5 * y[1]; dy)
        ts, ys = rkf45_integrate(rhs_exp!, [4.0], (0.0, 6.0); rtol = 1e-8, atol = 1e-10)
        @test isapprox(ys[1, end], 4.0 * exp(-3.0); rtol = 1e-6)
        @test ts[1] == 0.0 && isapprox(ts[end], 6.0; atol = 1e-12)

        # Harmonic oscillator energy conservation: E = 0.5 * (x^2 + v^2)
        rhs_osc!(dy, y, t) = (dy[1] = y[2]; dy[2] = -y[1]; dy)
        sample_t = collect(0.0:0.5:(4.0 * pi))
        ts_o, ys_o = rkf45_integrate(rhs_osc!, [1.0, 0.0], (0.0, 4.0 * pi); rtol = 1e-8, atol = 1e-10, sample_times = sample_t)
        @test length(ts_o) == length(sample_t)
        energies = [0.5 * (ys_o[1, k]^2 + ys_o[2, k]^2) for k in axes(ys_o, 2)]
        @test maximum(abs.(energies .- 0.5)) < 1e-5
    end

    @testset "Implicit Backward Euler & Crank-Nicolson Stiff Stability" begin
        # Stiff decay dy/dt = -50 * y
        rhs_stiff!(dy, y, t) = (dy[1] = -50.0 * y[1]; dy)
        _, y_ie = implicit_euler_integrate(rhs_stiff!, [2.0], (0.0, 1.0), 40)
        @test isfinite(y_ie[1, end]) && y_ie[1, end] >= 0.0 && y_ie[1, end] < 1e-4

        # Second-order convergence of Crank-Nicolson on dy/dt = -0.5 * y
        rhs_cn!(dy, y, t) = (dy[1] = -0.5 * y[1]; dy)
        _, y1 = crank_nicolson_integrate(rhs_cn!, [3.0], (0.0, 4.0), 20)
        _, y2 = crank_nicolson_integrate(rhs_cn!, [3.0], (0.0, 4.0), 40)
        exact = 3.0 * exp(-2.0)
        err1 = abs(y1[1, end] - exact)
        err2 = abs(y2[1, end] - exact)
        @test err1 / err2 > 3.7 # ~4x reduction for 2nd order
    end

    @testset "2D Poisson SOR Exact Manufactured Solution" begin
        nx, ny = 21, 21
        dx = 1.0 / (nx - 1)
        dy = 1.0 / (ny - 1)
        u = zeros(Float64, nx, ny)
        rhs = zeros(Float64, nx, ny)
        u_exact = zeros(Float64, nx, ny)
        for j in 1:ny, i in 1:nx
            x = (i - 1) * dx
            y = (j - 1) * dy
            u_exact[i, j] = sin(pi * x) * sin(pi * y)
            rhs[i, j] = -2.0 * pi^2 * u_exact[i, j]
        end
        iters, res = solve_poisson_sor!(u, rhs, dx, dy; tol = 1e-8, max_iter = 4000)
        @test res < 1e-7
        @test iters < 4000
        @test maximum(abs.(u .- u_exact)) < 3e-3
    end

    @testset "Brent Root, Levenberg-Marquardt, Welford & Morris" begin
        root = brent_root(x -> x^3 - 2.0 * x - 5.0, 1.0, 3.0; tol = 1e-11)
        @test isapprox(root^3 - 2.0 * root - 5.0, 0.0; atol = 1e-9)

        # Nonlinear least squares: y_i = p1 * exp(-p2 * t_i)
        ts = [0.0, 1.0, 2.0, 3.0, 4.0]
        y_obs = [3.5 * exp(-0.42 * t) for t in ts]
        p_opt, rmse, _, conv = levenberg_marquardt(
            p -> [p[1] * exp(-p[2] * ts[k]) - y_obs[k] for k in eachindex(ts)],
            [1.5, 0.1];
            lower = [0.1, 0.01],
            upper = [10.0, 2.0],
        )
        @test conv
        @test rmse < 1e-5
        @test isapprox(p_opt[1], 3.5; atol = 1e-4)
        @test isapprox(p_opt[2], 0.42; atol = 1e-4)

        acc = WelfordAccumulator()
        for v in [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
            push_sample!(acc, v)
        end
        @test isapprox(acc.mean, 5.0; atol = 1e-12)
        @test isapprox(sample_std(acc), std([2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]); atol = 1e-12)

        mu_star, _, _ = morris_sensitivity(x -> 10.0 * x[1] + 0.1 * x[2], [0.0, 0.0], [1.0, 1.0]; num_trajectories = 8)
        @test mu_star[1] > 20.0 * mu_star[2]
    end
end
