"""
Test integrating interpolated N-body dynamics with SPICE call within eom.
Uses low-level API
"""

using Base.Threads
using LinearAlgebra
using SPICE
using OrdinaryDiffEq
using Test

if !@isdefined(FullEphemerisPropagator)
    include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))
end

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
    f_jacobian = FullEphemerisPropagator.symbolic_Nbody_jacobian(length(mus)),
    naif_frame=naif_frame,
    abcorr=abcorr
)

et_range = (et0, et0 + 365.0*86400)
interp_params = FullEphemerisPropagator.InterpolatedNbody_params(et_range, parameters, 20000;
    rescale_epoch = false,)

# define state
u0_dim = [2200.0, 0.0, 42000.0, 0.03, 1.1, 0.1]
u0 = [u0_dim[1:3]/lstar; u0_dim[4:6]/parameters.vstar]

# propagation test with ensemble
tspan = (0.0, 4*86400/parameters.tstar)

# without STMs
N_traj = 5
x0_conditions = [
    [u0[1:3] + 10/parameters.lstar * randn(3); u0[4:6] + 1e-3/parameters.vstar * randn(3)]
    for _ in 1:N_traj
]
function prob_func_Nbody(ode_problem, i, repeat)
    _x0 = x0_conditions[i]
    remake(ode_problem, u0=_x0)
end

prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, parameters)
ensemble_prob = EnsembleProblem(
    prob;
    prob_func = prob_func_Nbody
)

prob_interp = ODEProblem(FullEphemerisPropagator.eom_Nbody_Interpolated!, u0, tspan, interp_params)
ensemble_prob_interp = EnsembleProblem(
    prob_interp;
    prob_func = prob_func_Nbody
)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial();
    trajectories=N_traj, reltol=1e-14, abstol=1e-14)

sols_interp = solve(ensemble_prob_interp, Vern9(), EnsembleThreads();
    trajectories=N_traj, reltol=1e-14, abstol=1e-14)

for (sol,sol_interp) in zip(sols,sols_interp)
    @test norm(sol.u[end] - sol_interp.u[end]) < 1e-11
    @test abs(sol.t[end] - sol_interp.t[end]) < 1e-16
end


# with STMs
N_traj = 5
x0_conditions = [
    [u0[1:3] + 10/parameters.lstar * randn(3); u0[4:6] + 1e-3/parameters.vstar * randn(3); reshape(I(6),36)]
    for _ in 1:N_traj
]
function prob_func_Nbody_STM(ode_problem, i, repeat)
    _x0 = x0_conditions[i]
    remake(ode_problem, u0=_x0)
end

prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_STM_SPICE!, u0, tspan, parameters)
ensemble_prob = EnsembleProblem(
    prob;
    prob_func = prob_func_Nbody_STM
)

prob_interp = ODEProblem(FullEphemerisPropagator.eom_Nbody_STM_Interpolated!, u0, tspan, interp_params)
ensemble_prob_interp = EnsembleProblem(
    prob_interp;
    prob_func = prob_func_Nbody_STM
)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial();
    trajectories=N_traj, reltol=1e-14, abstol=1e-14)

sols_interp = solve(ensemble_prob_interp, Vern9(), EnsembleThreads();
    trajectories=N_traj, reltol=1e-14, abstol=1e-14)

for (sol,sol_interp) in zip(sols,sols_interp)
    @test norm(sol.u[end][1:6] - sol_interp.u[end][1:6]) < 1e-11
    @test norm(sol.u[end][7:42] - sol_interp.u[end][7:42]) < 1e-10
    @test abs(sol.t[end] - sol_interp.t[end]) < 1e-16
end

