using Test
using LinearAlgebra
using Statistics
using IridisFlux

@testset "IridisFlux Full Test Suite" begin
    include("test_numerics.jl")
    include("test_aqueous.jl")
    include("test_anterior_pde.jl")
    include("test_biomechanics.jl")
    include("test_perfusion.jl")
    include("test_pkpd.jl")
    include("test_neuro_and_calibration.jl")
    include("test_cli_and_presets.jl")
end
