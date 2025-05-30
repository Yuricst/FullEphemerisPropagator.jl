"""Test ephemeris interpolation"""

using LinearAlgebra
using SPICE

if !@isdefined(FullEphemerisPropagator)
    include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))
end


test_interpolate_ephem = function ()
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


    interp_params = FullEphemerisPropagator.InterpolatedNbody_params(
        (et0, et0 + 30 * 86400.0),
        parameters,
        1000;
        rescale_epoch = false,
    )

    # evaluate position
    ets_test = et0 .+ LinRange(1e-3, 30 * 86400.0-1e-3, 1000)

    state_interp = zeros(6, length(ets_test))
    state_spice = zeros(6, length(ets_test))
    for (idx,et_test) in enumerate(ets_test)
        state_interp[:,idx] = FullEphemerisPropagator.get_state(interp_params.ephem_dict[2], et_test)

        state_spice[:,idx], _ = spkezr(
            naif_ids[2],
            et_test,
            naif_frame,
            abcorr,
            naif_ids[1]
        )
    end
    @test norm(state_spice[1:3,:] - state_interp[1:3,:] * lstar, Inf) < 1e-9
    @test norm(state_spice[4:6,:] - state_interp[4:6,:] * parameters.vstar, Inf) < 1e-12
    # @show norm(state_spice[1:3,:] - state_interp[1:3,:] * lstar, Inf)
    # @show norm(state_spice[4:6,:] - state_interp[4:6,:] * parameters.vstar, Inf)

    # # plot
    # components_str = ["x", "y", "z"]
    # diffs_pos = state_spice[1:3,:] - state_interp[1:3,:] * lstar
    # diffs_vel = (state_spice[4:6,:] - state_interp[4:6,:] * parameters.vstar)*1e6

    # fig = Figure(size=(1200, 600))
    # for i = 1:3
    #     ax_r = Axis(fig[1,i]; ylabel="Absolute error $(components_str[i]), km")
    #     lines!(ax_r, ets_test, diffs_pos[i,:])

    #     ax_v = Axis(fig[2,i]; ylabel="Absolute error v$(components_str[i]), mm/s")
    #     lines!(ax_v, ets_test, diffs_vel[i,:])
    # end
    # display(fig)
end

test_interpolate_ephem()