"""Parameters with interpolated ephemerides"""


abstract type InterpolatedEphemerisParameters end


"""
Construct parameters for interpolated equations of motion where no SPICE calls are made internally.
"""
mutable struct InterpolatedNbody_params <: InterpolatedEphemerisParameters
    et_range::Tuple{Float64, Float64}
    dt::Float64
    ephem_dict::Dict{Int, InterpolatedEphemeris}
    params::Nbody_params
    rescale_epoch::Bool

    et0::Float64
    lstar::Real
    tstar::Real
    vstar::Real
    mus_scaled::Vector{Float64}
    f_jacobian::Union{Nothing,Function}
    Rs

    function InterpolatedNbody_params(
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
        for i in 2:length(params.mus_scaled)
            # query states
            rvs_lt = spkezr.(
                params.naif_ids[i],
                ets,
                params.naif_frame,
                params.abcorr,
                params.naif_ids[1]
            )
            rvs = hcat([r_lt[1] for r_lt in rvs_lt]...)
            rvs[1:3,:] /= lstar_rescale
            rvs[4:6,:] /= vstar_rescale
             
            # fit interpolation
            ephem_dict[i] = InterpolatedEphemeris(
                params.naif_ids[i],
                ets,
                rvs,
                rescale_epoch,
                params.tstar
            )
        end
        new(et_range, ets[2] - ets[1], ephem_dict, params, rescale_epoch,
            params.et0, params.lstar, params.tstar, params.vstar, params.mus_scaled, params.f_jacobian,
            zeros(MVector{3 * (length(params.mus_scaled)-1),Float64}))
    end
end


"""
Overload method for showing InterpolatedNbody_params
"""
function Base.show(io::IO, itpparams::InterpolatedNbody_params)
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



"""
Construct parameters for interpolated equations of motion where no SPICE calls are made internally.
"""
mutable struct InterpolatedNbodySRP_params <: InterpolatedEphemerisParameters
    et_range::Tuple{Float64, Float64}
    dt::Float64
    ephem_dict::Dict{Int, InterpolatedEphemeris}
    params::NbodySRP_params
    rescale_epoch::Bool

    et0::Float64
    lstar::Real
    tstar::Real
    vstar::Real
    mus_scaled::Vector{Float64}
    k_srp::Real
    naif_ids::Vector{String}
    f_jacobian::Union{Nothing,Function}
    Rs

    function InterpolatedNbodySRP_params(
        et_range::Tuple{Float64, Float64},
        params::NbodySRP_params,
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
        for i in 2:length(params.mus_scaled)
            # query states
            rvs_lt = spkezr.(
                params.naif_ids[i],
                ets,
                params.naif_frame,
                params.abcorr,
                params.naif_ids[1]
            )
            rvs = hcat([r_lt[1] for r_lt in rvs_lt]...)
            rvs[1:3,:] /= lstar_rescale
            rvs[4:6,:] /= vstar_rescale
             
            # fit interpolation
            ephem_dict[i] = InterpolatedEphemeris(
                params.naif_ids[i],
                ets,
                rvs,
                rescale_epoch,
                params.tstar
            )
        end
        new(et_range, ets[2] - ets[1], ephem_dict, params, rescale_epoch,
            params.et0, params.lstar, params.tstar, params.vstar,
            params.mus_scaled, params.k_srp, params.naif_ids, params.f_jacobian,
            zeros(MVector{3 * (length(params.mus_scaled)-1),Float64}))
    end
end