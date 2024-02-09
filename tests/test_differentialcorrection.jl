"""
Test jacobian function
"""

using SPICE
using GLMakie
using DifferentialEquations

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

# furnish spice kernels
spice_dir = ENV["SPICE"]

# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
furnsh(joinpath(spice_dir, "fk", "earth_moon_rotating_mc.tf"))

naif_frame = "J2000"
lstar = 3000
mus = [
    4.9028000661637961E+03,
    3.9860043543609598E+05,
    1.3271244004193938E+11,
]
naif_ids = ["301", "399", "10"]
N = length(naif_ids)

# initialize integrator
prop = FullEphemerisPropagator.PropagatorSTM(
    Vern9(),
    lstar,
    mus,
    naif_ids;
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# initial epoch
et0 = str2et("2030-04-16T00:00:00")

# initial guess in CR3BP
LU = 389703.0
TU = 382981.0
mu_cr3bp = 1.215058560962404E-2
# DRO ID 943
x0_cr3bp = [
    9.4916844627421393E-1, 0.0, 0.0, 
    0.0, 6.0205794060789086E-1, 0.0
] - [1-mu_cr3bp, 0.0, 0.0, 0.0, 0.0, 0.0]   # shifted to be Moon centered
period_cr3bp = 4.0938520066556927E-1

x0_EMrot = vcat(x0_cr3bp[1:3]*LU, x0_cr3bp[4:6]*LU/TU)
T_EM2Inr = sxform("EARTHMOONROTATINGMC", naif_frame, et0)
x0_SI = T_EM2Inr * x0_EMrot

# convert to scale of integrator
#x0 = FullEphemerisPropagator.dim2nondim(prop, x0_SI)
x0 = vcat(x0_SI[1:3]/prop.parameters.lstar, x0_SI[4:6]/prop.parameters.vstar)

# propagate
tspan = (0.0, period_cr3bp * TU / prop.parameters.tstar)
sol = FullEphemerisPropagator.propagate(prop, et0, x0, tspan)

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
lines!(ax1, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:])
nsph = 30
θ = range(0, stop=2π, length=nsph)
ϕ = range(0, stop=π, length=nsph)
R = 1737.4/lstar
center = [0, 0, 0]
xsphere = [center[1] + R * cos(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
ysphere = [center[2] + R * sin(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
zsphere = [center[3] + R * cos(ϕ[j]) for j in 1:nsph, i in 1:nsph]
wireframe!(ax1, xsphere, ysphere, zsphere, color=:grey, linewidth=0.5)


# construct differential correction problem
Nrev = 3
epochs = [et0 + (idx-1) * period_cr3bp * TU for idx in 1:Nrev]
nodes = [x0 for _ in 1:Nrev]

problem = FullEphemerisPropagator.ForwardMultipleShootingProblem(
    prop,
    epochs,
    nodes,
)

# solve
sols, residuals, DF = FullEphemerisPropagator.shoot(problem, maxiter=2);

for _sol in sols
    lines!(ax1, Array(_sol)[1,:], Array(_sol)[2,:], Array(_sol)[3,:])
end


display(fig)
println("Done!")


