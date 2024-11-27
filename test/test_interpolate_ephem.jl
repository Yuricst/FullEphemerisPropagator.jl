"""Test ephemeris interpolation"""



include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))


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


interp_params = FullEphemerisPropagator.InterpolatedNbodyParams(
    (et0, et0 + 30 * 86400.0),
    parameters,
    1000;
    rescale_epoch = false,
)

# evaluate position
ets_test = et0 .+ LinRange(1e-3, 30 * 86400.0-1e-3, 1000)

pos_interp = zeros(3, length(ets_test))
pos_spice = zeros(3, length(ets_test))
for (idx,et_test) in enumerate(ets_test)
    pos_interp[:,idx] = FullEphemerisPropagator.get_pos(interp_params.ephem_dict[2], et_test) * lstar

    pos_spice[:,idx], _ = spkpos(
        naif_ids[2],
        et_test,
        naif_frame,
        abcorr,
        naif_ids[1]
    )
end
@show norm(pos_spice - pos_interp, Inf)

# plot
components_str = ["x", "y", "z"]
diffs = pos_spice - pos_interp
fig = Figure(size=(1200, 600))
for i = 1:3
    ax_r = Axis(fig[1,i]; ylabel="Absolute error $(components_str[i]), km")
    lines!(ax_r, ets_test, diffs[i,:])
end
display(fig)