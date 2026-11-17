@testset "CLI Subcommands, Presets & Sensitivity Analysis" begin
    @testset "All Clinical Phenotypes & Drug Regimens Valid" begin
        for pname in ("healthy", "poag", "ntg", "aacc", "pxfg")
            p = get_phenotype(pname)
            @test p.anatomy.globe_radius_mm > 0.0
            @test steady_state_iop(p.aqueous; hour_of_day = 10.0) > 8.0
        end
        for rname in ("none", "latanoprost", "timolol", "dorzolamide", "brimonidine", "netarsudil", "dual", "max")
            r = get_regimen(rname)
            @test r.adherence_probability > 0.0
        end
        @test_throws ArgumentError get_phenotype("nonexistent")
        @test_throws ArgumentError get_regimen("nonexistent")
    end

    @testset "Morris Sensitivity Ranking" begin
        poag = poag_moderate_phenotype()
        sens = run_aqueous_sensitivity(poag; num_trajectories = 6)
        @test length(sens) == 5
        @test issorted([r.mu_star for r in sens]; rev = true)
        @test sens[1].mu_star > 0.0
    end

    @testset "CLI Dispatcher Subcommands" begin
        for argv in (
            ["--help"],
            ["simulate", "--phenotype", "poag", "--drug", "latanoprost", "--hours", "12", "--plot"],
            ["anterior", "--steps", "30", "--plot"],
            ["calibrate", "--phenotype", "poag"],
            ["benchmark"],
            ["inspect", "--phenotype", "ntg"],
        )
            buf = IOBuffer()
            code = run_cli(String.(argv); io = buf)
            out = String(take!(buf))
            @test code == 0
            @test !isempty(out)
            @test all(c -> UInt32(c) <= 127, out)
        end
        buf_err = IOBuffer()
        @test run_cli(["unknown_cmd"]; io = buf_err) == 1
    end
end
