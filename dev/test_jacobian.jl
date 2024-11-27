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

naif_frame = "J2000"
lstar = 3000
mus = [
    4.9028000661637961E+03,
    3.9860043543609598E+05,
    1.3271244004193938E+11,
]
naif_ids = ["301", "399", "10"]
N = length(naif_ids)

# # test function
# println("Generating Jacobian function symbolically...")
# @time f_jacobian = FullEphemerisPropagator.symbolic_Nbody_jacobian(N);
# Rs = hcat(SPICE.spkpos("399", 0.0, "J2000", "NONE", "301")[1], 
#           SPICE.spkpos("10", 0.0, "J2000", "NONE", "301")[1])
# f_jacobian([u0_dim[1:3]..., mus..., Rs])

# integrate with STM
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
et0 = 840463475.5117936 # str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
# u0_dim = [3200.0, 0.0, 4200.0, 0.0, 0.77, 0.0]
# u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)
u0 = [
    -2.5019204591096096,
    14.709398066624694,
    -18.59744250295792,
    5.62688812721852e-2,
    1.439926311669468e-2,
    3.808273517470642e-3
]

# time span (in canonical scale)
tspan = (0.0, 6.55 * 5 *86400/prop.parameters.tstar)

# evaluate jacobian at initial time
jac = FullEphemerisPropagator.jacobian(prop, et0, 0.0, u0)
println("jac:")
display(jac)

# solve
println("Integrating...")
@time sol = FullEphemerisPropagator.propagate(prop, et0, tspan, u0)
println("Done!")
@show sol.u[end][1:6], sol.t[end]
display(reshape(sol.u[end][7:42], (6,6)))

# validate STM
prop_no_stm = FullEphemerisPropagator.Propagator(
    Vern9(),
    lstar,
    mus,
    naif_ids;
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)
STM_numerical = zeros(6,6)
for i = 1:6
    u0_ptrb = deepcopy(u0)
    u0_ptrb[i] += 1e-6
    _sol = FullEphemerisPropagator.propagate(prop_no_stm, et0, tspan, u0_ptrb)
    STM_numerical[:,i] = (_sol.u[end][1:6] - sol.u[end][1:6]) / 1e-6
end
display(STM_numerical)

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