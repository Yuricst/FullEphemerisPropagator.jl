"""
Test jacobian function
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
furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))
furnsh(joinpath(spice_dir, "fk", "earth_moon_rotating_mc.tf"))

# define parameters
naif_ids = ["301", "399", "10"]
mus = [bodvrd(ID, "GM", 1)[1] for ID in naif_ids]
naif_frame = "J2000"
abcorr = "NONE"
lstar = 3000.0

# initialize integrator
prop = FullEphemerisPropagator.Propagator(
    Vern9(),
    lstar,
    mus,
    naif_ids;
    use_srp = true,
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)
prop_stm = FullEphemerisPropagator.PropagatorSTM(
    Vern9(),
    lstar,
    mus,
    naif_ids;
    use_srp = true,
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

# create initial guess
Nrev = 10
epochs = Float64[]
nodes = Vector{Float64}[]
for irev in 1:Nrev
    # compue state at epoch, in CR3BP canonical scales
    et_rev = et0 + (irev-1) * period_cr3bp * TU
    T_EM2Inr = sxform("EARTHMOONROTATINGMC", naif_frame, et_rev)
    x0_SI = T_EM2Inr * x0_EMrot

    # convert to scale of integrator
    node = vcat(x0_SI[1:3]/prop_stm.parameters.lstar, x0_SI[4:6]/prop_stm.parameters.vstar)

    push!(epochs, et_rev)
    push!(nodes, node)
end

# construct differential correction problem
problem = FullEphemerisPropagator.ForwardMultipleShootingProblem(
    prop,
    prop_stm,
    epochs,
    nodes,
)
@show problem.epochs
@show problem.tofs

# solve
sols, residuals, DF = FullEphemerisPropagator.shoot_freetime(problem, maxiter=5);
@show problem.epochs
@show problem.tofs

# plot with GLMakie
fig = Figure(size=(600,400), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=:data)
nsph = 30
θ = range(0, stop=2π, length=nsph)
ϕ = range(0, stop=π, length=nsph)
R = 1737.4/lstar
center = [0, 0, 0]
xsphere = [center[1] + R * cos(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
ysphere = [center[2] + R * sin(θ[i]) * sin(ϕ[j]) for j in 1:nsph, i in 1:nsph]
zsphere = [center[3] + R * cos(ϕ[j]) for j in 1:nsph, i in 1:nsph]
wireframe!(ax1, xsphere, ysphere, zsphere, color=:grey, linewidth=0.5)

for _sol in sols
    lines!(ax1, Array(_sol)[1,:], Array(_sol)[2,:], Array(_sol)[3,:])
    scatter!(ax1, Array(_sol)[1,1], Array(_sol)[2,1], Array(_sol)[3,1], color=:green)
    scatter!(ax1, Array(_sol)[1,end], Array(_sol)[2,end], Array(_sol)[3,end], color=:red, marker=:rect)
end

display(fig)
println("Done!")


