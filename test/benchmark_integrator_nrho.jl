"""Run benchmark on integrator"""

using SPICE
using GLMakie
using OrdinaryDiffEq
using StaticArrays
using BenchmarkTools
using ProgressMeter: @showprogress
using Printf: @printf

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

# furnish spice kernels
spice_dir = ENV["SPICE"]

# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))

# define parameters
naif_ids = ["301", "399", "10"]
mus = [bodvrd(ID, "GM", 1)[1] for ID in naif_ids]
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

prop_stm = FullEphemerisPropagator.PropagatorSTM(
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
et0 = 788991045.1839374
if use_sa
    u0 = SA[5.0883176857384, 0.39296133080885437, -19.251718764456303,
        -0.03196443660603664, -0.015104367539665122, -0.15169141145203588]
else
    u0 = [5.0883176857384, 0.39296133080885437, -19.251718764456303,
        -0.03196443660603664, -0.015104367539665122, -0.15169141145203588]
end

# time span (in canonical scale)
tf_day = 30.0
tspan = (0.0, tf_day*86400/prop.parameters.tstar)

# solve state only
N_sample = 10
t_measure = []
@showprogress for i = 1:N_sample
    tstart = time()
    FullEphemerisPropagator.propagate(prop, et0, tspan, u0);
    push!(t_measure, time() - tstart)
end
println("FullEphemerisPropagator.jl Propagator over $(tf_day) days ($(N_sample) run samples)")
@printf("    mean exec. time : %1.4f s\n", mean(t_measure))
if N_sample > 1
    @printf("    std dev. time   : %1.4f s\n", std(t_measure))
end

# solve state and STM
N_sample = 10
t_measure = []
@showprogress for i = 1:N_sample
    tstart = time()
    FullEphemerisPropagator.propagate(prop_stm, et0, tspan, u0);
    push!(t_measure, time() - tstart)
end
println("FullEphemerisPropagator.jl PropagatorSTM over $(tf_day) days ($(N_sample) run samples)")
@printf("    mean exec. time : %1.4f s\n", mean(t_measure))
if N_sample > 1
    @printf("    std dev. time   : %1.4f s\n", std(t_measure))
end

println("Done!")
# println("Benchmarking...")
# @benchmark sol = FullEphemerisPropagator.propagate(prop, et0, tspan, u0; saveat=tspan[2])
