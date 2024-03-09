"""
Test integrating N-body dynamics with SPICE call within eom
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

# instantiate propagator
prop = FullEphemerisPropagator.Propagator(
    Vern7(),
    lstar,
    mus,
    naif_ids;
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
u0_dim = [3200.0, 0.0, 4200.0, 0.0, 0.77, 0.0]
u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)

# time span (in canonical scale)
tspan = (0.0, 10*86400/prop.parameters.tstar)

# solve
sol = FullEphemerisPropagator.propagate(prop, et0, u0, tspan)
@show sol.u[end];

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])

# Generate points on the sphere
nsph = 30
θ = range(0, stop=2π, length=nsph)
ϕ = range(0, stop=π, length=nsph)
R = 1737.4/lstar
center = [0, 0, 0]
xsphere = [center[1] + R * cos(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
ysphere = [center[2] + R * sin(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
zsphere = [center[3] + R * cos(ϕ[j]) for j in 1:nsph, i in 1:nsph]
wireframe!(ax1, xsphere, ysphere, zsphere, color=:grey, linewidth=0.5)

fig