"""Run benchmark on integrator"""

using SPICE
using GLMakie
using OrdinaryDiffEq
using StaticArrays
using BenchmarkTools

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

# furnish spice kernels
spice_dir = ENV["SPICE"]

# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# define parameters
mus = [
    4.9028000661637961E+03,
    3.9860043543609598E+05,
    1.3271244004193938E+11,
]
naif_ids = ["301", "399", "10"]
naif_frame = "J2000"
abcorr = "NONE"
lstar = 3000.0

# instantiate propagator
use_sa = true
prop = FullEphemerisPropagator.Propagator(
    Vern7(),
    lstar,
    mus,
    naif_ids;
    use_srp = false,
    use_sa = use_sa,
    naif_frame = naif_frame,
    reltol = 1e-14,
    abstol = 1e-12,
)

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
if use_sa == true
    u0_dim = SA[2200.0, 0.0, 4200.0, 0.03, 1.1, 0.1]
else
    u0_dim = [2200.0, 0.0, 4200.0, 0.03, 1.1, 0.1]
end
u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)

# time span (in canonical scale)
tspan = (0.0, 1*86400/prop.parameters.tstar)

# solve
sol = FullEphemerisPropagator.propagate(prop, et0, tspan, u0)
@show sol.u[end][1:6]

# println("Benchmarking...")
# @benchmark sol = FullEphemerisPropagator.propagate(prop, et0, tspan, u0; saveat=tspan[2])
