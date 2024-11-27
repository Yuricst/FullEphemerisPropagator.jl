"""Test STM propagation"""

using SPICE
using GLMakie
using OrdinaryDiffEq
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

# instantiate propagator without STM
prop = FullEphemerisPropagator.Propagator(
    Vern7(),
    lstar,
    mus,
    naif_ids;
    use_srp = true,
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# instantiate propagator with STM
prop_stm = FullEphemerisPropagator.PropagatorSTM(
    Vern7(),
    lstar,
    mus,
    naif_ids;
    use_srp = true,
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
u0_dim = [2200.0, 0.0, 4200.0, 0.03, 1.1, 0.1]
u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)

# time span (in canonical scale)
tspan = (0.0, 1*86400/prop.parameters.tstar)

# solve with STM
sol_stm = FullEphemerisPropagator.propagate(prop_stm, et0, tspan, u0)
@show sol_stm.u[end][1:6]
println("STM via propagation:")
display(reshape(sol_stm.u[end][7:end], 6, 6)')

# solve state only & get STM via forward differencing
h_step = 1e-6
stm_num = zeros(6,6)
sol_nominal = FullEphemerisPropagator.propagate(prop, et0, tspan, u0)
@show sol_nominal.u[end][1:6]
for i = 1:6
    _u0 = copy(u0)
    _u0[i] += h_step
    sol_ptrb = FullEphemerisPropagator.propagate(prop, et0, tspan, _u0)
    stm_num[:,i] = (sol_ptrb.u[end] - sol_nominal.u[end]) / h_step
end
println("STM via forward difference:")
display(stm_num)