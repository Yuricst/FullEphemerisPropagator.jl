"""Interpolate ephemeris"""


struct InterpolatedEphemeris
    naif_id::String
    et_range::Tuple{Float64, Float64}
    splines::Array{Spline1D, 1}
    rescale_epoch::Bool
    tstar::Float64

    function InterpolatedEphemeris(
        naif_id::String,
        ets,
        rvs,
        rescale_epoch::Bool,
        tstar::Float64
    )
        if rescale_epoch
            @warn "rescale_epoch == true is buggy"
            times_input = (ets .- ets[1]) / tstar
        else
            times_input = ets
        end
        splines = [
            Spline1D(times_input, rvs[1,:]),
            Spline1D(times_input, rvs[2,:]),
            Spline1D(times_input, rvs[3,:]),
            Spline1D(times_input, rvs[4,:]),
            Spline1D(times_input, rvs[5,:]),
            Spline1D(times_input, rvs[6,:]),
        ]
        new(naif_id, (ets[1], ets[end]), splines, rescale_epoch, tstar)
    end
end


"""
Overload method for showing InterpolatedEphemeris
"""
function Base.show(io::IO, ephem::InterpolatedEphemeris)
    println("Interpolated ephemeris struct")
    @printf("    et0        : %s (%1.8f)\n", et2utc(ephem.et_range[1], "ISOC", 3), ephem.et_range[1])
    @printf("    etf        : %s (%1.8f)\n", et2utc(ephem.et_range[2], "ISOC", 3), ephem.et_range[2])
    @printf("    naif_id    : %s\n", ephem.naif_id)
end


"""Interpolate ephemeris position at a given epoch"""
function get_pos(ephem::InterpolatedEphemeris, et::Float64)
    if ephem.rescale_epoch
        et_eval = et * ephem.tstar + ephem.et_range[1]
        @assert ephem.et_range[1] <= et <= ephem.et_range[2]
    else
        et_eval = et
        @assert ephem.et_range[1] <= et <= ephem.et_range[2]
    end
    return [Dierckx.evaluate(ephem.splines[1], et_eval),
            Dierckx.evaluate(ephem.splines[2], et_eval),
            Dierckx.evaluate(ephem.splines[3], et_eval)]
end


"""Interpolate ephemeris state at a given epoch"""
function get_state(ephem::InterpolatedEphemeris, et::Float64)
    if ephem.rescale_epoch
        et_eval = et * ephem.tstar + ephem.et_range[1]
        @assert ephem.et_range[1] <= et <= ephem.et_range[2]
    else
        et_eval = et
        @assert ephem.et_range[1] <= et <= ephem.et_range[2]
    end
    return [Dierckx.evaluate(ephem.splines[1], et_eval),
            Dierckx.evaluate(ephem.splines[2], et_eval),
            Dierckx.evaluate(ephem.splines[3], et_eval),
            Dierckx.evaluate(ephem.splines[4], et_eval),
            Dierckx.evaluate(ephem.splines[5], et_eval),
            Dierckx.evaluate(ephem.splines[6], et_eval)]
end

