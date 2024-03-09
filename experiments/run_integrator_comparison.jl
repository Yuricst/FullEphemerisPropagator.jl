"""
Run integrator benchmarking comparison
"""

using GLMakie
using SPICE
using OrdinaryDiffEq
using BenchmarkTools
using ProgressMeter
using LinearAlgebra

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))

function forwardbackward(et0, tf, x0)#, steps)
    # 1. propagate forward
    sol_fwd = FullEphemerisPropagator.propagate(prop, et0, (0.0, tf), x0
        #saveat=LinRange(0.0, tf, steps)
    )
    # 2. propagate backward
    etf = et0 + FullEphemerisPropagator.TU2sec(prop, tf)
    sol_bck = FullEphemerisPropagator.propagate(
        prop, etf, (0.0, -tf), sol_fwd.u[end]
        #saveat=LinRange(0.0, -tf, steps)
    )
    # compute error
    error = sol_bck.u[end] - x0
    return error, sol_fwd, sol_bck
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

# instantiate propagator
prop = FullEphemerisPropagator.Propagator(
    Vern8(),
    lstar,
    mus,
    naif_ids;
    use_srp = true,
    naif_frame = naif_frame,
    reltol = 1e-14,
    abstol = 1e-12,
);

# ---------------------------------------------------------------------- #
orbit_type = "NRHO"  # NRHO or LO15000
@show orbit_type

# initial epoch and state of an NRHO (almost) in canonical scales using LU = 3000 km
et0 = 946728069.3271508
if orbit_type == "NRHO"
    state0 = [-0.03344377115230989, 5.7624151996473545, -22.743438743043676, 
              -0.046526421787245704, 0.029709480647552224, 0.004309142532513644]
elseif orbit_type == "LO15000"
    R_LLO = 1737.4 + 15000.0
    state0 = [R_LLO/lstar, 0.0, 0.0,
              0.0, 0.0, sqrt(mus[1]/R_LLO)/prop.parameters.vstar]
end

#steps = 500
tspan = (0.0, FullEphemerisPropagator.sec2TU(prop, 30 * 86400))
#tevals = LinRange(tspan[1], tspan[2], steps)

# methods to be used
methods = [Vern7(), TanYam7(),            # RK 7
           Vern8(), TsitPap8(), DP8(),    # RK 8(7) methods
           Vern9(),                       # RK 9(8) methods                  
]
FullEphemerisPropagator.propagate(prop, et0, (0.0, 1.0), state0)  # run for jit-compiling

# ---------------------------------------------------------------------- #
# run speed test
println("********** 1. Speed test **********")
bruns = []
@showprogress for method in methods
    # mutate method of propagator
    prop.method = method
    
    # dry run once to ensure it is pre-compiled
    FullEphemerisPropagator.propagate(prop, et0, tspan, state0);
    
    # run benchmark
    #FullEphemerisPropagator.pretty(prop)
    b = @benchmarkable FullEphemerisPropagator.propagate($prop,
        $et0, $tspan, $state0; #saveat=$tevals
    );
    tune!(b);
    push!(bruns, run(b))
end

# ---------------------------------------------------------------------- #
# run accuracy test
println("********** 2. Accuracy test **********")
tf = FullEphemerisPropagator.sec2TU(prop, 30 * 86400)
errors = []
n_timesteps = Int[]

@showprogress for method in methods
    # mutate method of propagator
    prop.method = method

    # run forward-backward test
    error, sol_fwd, _ = forwardbackward(et0, tf, state0, steps)
    push!(errors, [norm(error[1:3]), norm(error[4:6])])
    push!(n_timesteps, length(sol_fwd.t))
end
errors = hcat(errors...)

# ---------------------------------------------------------------------- #
# make plot
fontsize = 18
names = [split(string(method), "(;")[1] for method in methods]
fig = Figure(size=(1000,800))

xlabelrot = 60

# make histogram of times
ax = Axis(fig[1:3,1];
    ylabel = "Integration time, ms", xticks = (1:length(names), names),
    xticklabelrotation = deg2rad(xlabelrot),
    titlesize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize-1, yticklabelsize=fontsize-1,
)
categories = vcat([[idx for _ in 1:length(brun.times)] for (idx,brun) in enumerate(bruns)]...,)
values = vcat([brun.times/1e9 for brun in bruns]...)
violin!(ax, categories, values, color=:navy)

# make scatters of number of time-steps taken
ax = Axis(fig[1,2];
    ylabel = "Time-steps", xticks = (1:length(names), names), #yscale = log10,
    xticklabelrotation = deg2rad(xlabelrot),
    titlesize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize-1, yticklabelsize=fontsize-1,
)
scatter!(ax, 1:length(names), n_timesteps, color=:black)

# make scatter of errors
ax = Axis(fig[2,2];
    ylabel = "Position error, m", xticks = (1:length(names), names), yscale = log10,
    xticklabelrotation = deg2rad(xlabelrot),
    titlesize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize-1, yticklabelsize=fontsize-1,
)
scatter!(ax, 1:length(names), errors[1,:] * prop.parameters.lstar * 1e3, color=:crimson)
ax = Axis(fig[3,2];
    ylabel = "Velocity error, mm/s", xticks = (1:length(names), names), yscale = log10,
    xticklabelrotation = deg2rad(xlabelrot),
    titlesize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize-1, yticklabelsize=fontsize-1,
)
scatter!(ax, 1:length(names), errors[2,:] * prop.parameters.vstar * 1e6, color=:crimson)

supertitle = Label(fig[0,:], "Integration of $(orbit_type) over 30 days (reltol = 1e-14, abstol = 1e-12)", fontsize = fontsize+2)
save(joinpath(dirname(@__FILE__), "integrator_comparison_$(orbit_type).png"), fig)
fig
