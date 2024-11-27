"""Run tests"""

using LinearAlgebra
using OrdinaryDiffEq
using SPICE
using Test

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

# furnish spice kernels
spice_dir = ENV["SPICE"]
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))

@time @testset "SPICE-based Propagators" begin
    include("test_jacobian.jl")
    include("test_Nbody_spice.jl")
    include("test_NbodySRP_spice.jl")
    include("test_Nbody_spice_lowapi.jl")
    include("test_stm.jl")
end

@time @testset "Ephemeris interpolation" begin
    include("test_interpolate_ephem.jl")
end