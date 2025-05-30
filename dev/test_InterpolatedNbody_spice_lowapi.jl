"""
Test integrating interpolated N-body dynamics with SPICE call within eom.
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
    f_jacobian = FullEphemerisPropagator.symbolic_Nbody_jacobian(length(mus)),
    naif_frame=naif_frame,
    abcorr=abcorr
)

et_range = (et0, et0 + 365.0*86400)
interp_params = FullEphemerisPropagator.InterpolatedNbody_params(et_range, parameters, 20000;
     rescale_epoch = false,)

# define state
u0_dim = [2200.0, 0.0, 4200.0, 0.03, 1.1, 0.1]
u0 = [u0_dim[1:3]/lstar; u0_dim[4:6]/parameters.vstar]

drv = FullEphemerisPropagator.eom_Nbody_SPICE(u0, parameters, 0.0)
drv_interp = FullEphemerisPropagator.eom_Nbody_SPICE(u0, interp_params, 0.0)
@test norm(drv - drv_interp) < 1e-14

drv_inplace_interp = zeros(6)
FullEphemerisPropagator.eom_Nbody_SPICE!(drv_inplace_interp, u0, interp_params, 0.0)
@test norm(drv - drv_inplace_interp) < 1e-14

u0_stm = [u0; reshape(I(6),36)]
drvstm_inplace = zeros(42)
drvstm_inplace_interp = zeros(42)
FullEphemerisPropagator.eom_Nbody_STM_SPICE!(drvstm_inplace, u0_stm, interp_params, 0.0)
FullEphemerisPropagator.eom_Nbody_STM_SPICE!(drvstm_inplace_interp, u0_stm, interp_params, 0.0)
@test norm(drvstm_inplace - drvstm_inplace_interp) < 1e-14

# propagation test
tspan = (0.0, 7*86400/parameters.tstar)

# solve with SPICE-based parameters
tstart = time()
prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, parameters)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
println("SPICE-based propagation time: $(time() - tstart) seconds")
@show sol.u[end];

# solve with interpolated parameters
tstart = time()
prob_interp = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, interp_params)
sol_interp = solve(prob_interp, Tsit5(), reltol=1e-12, abstol=1e-12)
println("Interpolated propagation time: $(time() - tstart) seconds")
@show sol_interp.u[end];

errors = Array(sol_interp) - hcat(sol.(sol_interp.t)...)
@test norm(errors) < 1e-10

# plot with GLMakie
fig = Figure(size=(900,900))
ax1 = Axis3(fig[1, 1], aspect=(1,1,1), xlabel="x", ylabel="y", zlabel="z", title="SPICE-based")
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])

ax2 = Axis3(fig[1, 2], aspect=(1,1,1), xlabel="δx", ylabel="δy", zlabel="δz", title="Interpolated")
lines!(ax2, sol_interp[1,:], sol_interp[2,:], sol_interp[3,:])

ax3 = Axis3(fig[2,1], aspect=(1,1,1), xlabel="δx", ylabel="δy", zlabel="δz")
lines!(ax3, errors[1,:], errors[2,:], errors[3,:])
scatter!(ax3, errors[1,end], errors[2,end], errors[3,end], marker=:circle)

ax4 = Axis3(fig[2,2], aspect=(1,1,1), xlabel="δvx", ylabel="δvy", zlabel="δvz")
lines!(ax4, errors[4,:], errors[5,:], errors[6,:])
scatter!(ax4, errors[4,end], errors[5,end], errors[6,end], marker=:circle)
fig