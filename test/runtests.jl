"""Run tests"""

using GLMakie
using LinearAlgebra
using OrdinaryDiffEq
using SPICE
using Test

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

# furnish spice kernels
spice_dir = ENV["SPICE"]
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))



@time @testset "Propagator object tests" begin
    include("test_stm.jl")
end