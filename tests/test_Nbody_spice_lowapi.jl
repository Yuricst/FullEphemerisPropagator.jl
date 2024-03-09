"""
Test integrating N-body dynamics with SPICE call within eom.
Uses low-level API
"""

using SPICE
using GLMakie
using OrdinaryDiffEq

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

et0 = str2et("2020-01-01T00:00:00")
parameters = FullEphemerisPropagator.Nbody_params(
    et0,
    lstar,
    mus,
    naif_ids;
    naif_frame=naif_frame,
    abcorr=abcorr
)

# initial state (in canonical scale)
u0 = [1.0, 0.0, 0.3, 0.0, 1.0, 0.0]

# time span (in canonical scale)
tspan = (0.0, 7*86400/parameters.tstar)

# solve
prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, parameters)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
@show sol.u[end];

# plot with GLMakie
fig = Figure(resolution=(600,600), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=(1,1,1))
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])
fig