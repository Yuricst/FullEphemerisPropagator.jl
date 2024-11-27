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


abstract type InterpolatedEphemerisParameters end

mutable struct InterpolatedNbodyParams <: InterpolatedEphemerisParameters
    et_range::Tuple{Float64, Float64}
    dt::Float64
    ephem_dict::Dict{Int, InterpolatedEphemeris}
    params::Nbody_params
    rescale_epoch::Bool

    function InterpolatedNbodyParams(
        et_range::Tuple{Float64, Float64},
        params::Nbody_params,
        N_step::Int = 10;
        rescale_epoch::Bool = false,
        rescale_state::Bool = true,
    )
        if rescale_epoch
            @warn "rescale_epoch == true is buggy"
        end

        # scaling
        # scaling
        if rescale_state
            lstar_rescale = params.lstar
            vstar_rescale = params.vstar
        else
            lstar_rescale = 1.0
            vstar_rescale = 1.0
        end

        # create splines
        ephem_dict = Dict{Int, InterpolatedEphemeris}()
        ets = range(et_range[1], et_range[2], N_step)
        for i = 2:length(params.mus_scaled)
            rvs_lt = spkezr.(
                params.naif_ids[2],
                ets,
                params.naif_frame,
                params.abcorr,
                params.naif_ids[1]
            )
            rvs = hcat([r_lt[1] for r_lt in rvs_lt]...)
            rvs[1:3,:] /= lstar_rescale
            rvs[4:6,:] /= vstar_rescale

            ephem_dict[i] = InterpolatedEphemeris(
                params.naif_ids[2],
                ets,
                rvs,
                rescale_epoch,
                params.tstar
            )
        end
        new(et_range, ets[2] - ets[1], ephem_dict, params, rescale_epoch)
    end
end


"""
Overload method for showing InterpolatedNbodyParams
"""
function Base.show(io::IO, itpparams::InterpolatedNbodyParams)
    println("Interpolated ephemeris parameter struct")
    @printf("    et0              : %s (%1.8f)\n", et2utc(itpparams.et_range[1], "ISOC", 3), itpparams.et_range[1])
    @printf("    etf              : %s (%1.8f)\n", et2utc(itpparams.et_range[2], "ISOC", 3), itpparams.et_range[2]) 
    @printf("    interpolation dt : %1.4f sec\n", itpparams.dt)
    @printf("    lstar            : %1.4f\n", itpparams.params.lstar)
    @printf("    tstar            : %1.4f\n", itpparams.params.tstar)
    @printf("    vstar            : %1.4f\n", itpparams.params.vstar)
    @printf("    mus_scaled       : %s\n", string(itpparams.params.mus_scaled))
    @printf("    naif_ids         : %s\n", string(itpparams.params.naif_ids))
    @printf("    naif_frame       : %s\n", itpparams.params.naif_frame)
    @printf("    abcorr           : %s\n", itpparams.params.abcorr)
end