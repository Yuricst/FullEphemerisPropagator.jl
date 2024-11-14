"""
Parameters file
"""

abstract type FullEphemParameters end

mutable struct Nbody_params <: FullEphemParameters
    et0::Float64
    lstar::Real
    tstar::Real
    vstar::Real
    mus::Vector{Float64}
    mus_scaled::Vector{Float64}
    naif_ids::Vector{String}
    naif_frame::String
    abcorr::String
    f_jacobian::Union{Nothing,Function}
    Rs

    # constructor
    function Nbody_params(
        et0,
        lstar::Real,
        mus::Vector{Float64},
        naif_ids::Vector{String};
        naif_frame::String="J2000",
        abcorr::String="NONE",
        use_canonical::Bool = true,
    )
        if use_canonical
            # scaled mus
            mus_scaled = mus / mus[1]

            # calculate tstar and vstar
            vstar = sqrt(mus[1]/lstar)
            tstar = lstar/vstar
        else
            mus_scaled = mus
            lstar = 1.0
            vstar = 1.0
            tstar = 1.0
        end

        # initialize object
        new(
            et0,
            lstar,
            tstar,
            vstar,
            mus,
            mus_scaled,
            naif_ids,
            naif_frame,
            abcorr,
            nothing,
            zeros(MVector{3 * (length(mus)-1),Float64})
        )
    end
end


mutable struct NbodySRP_params <: FullEphemParameters
    et0::Float64
    lstar::Real
    tstar::Real
    vstar::Real
    mus::Vector{Float64}
    mus_scaled::Vector{Float64}
    k_srp::Real
    naif_ids::Vector{String}
    naif_frame::String
    abcorr::String
    f_jacobian::Union{Nothing,Function}
    Rs

    # constructor
    function NbodySRP_params(
        et0,
        lstar::Real,
        mus::Vector{Float64},
        naif_ids::Vector{String},
        srp_cr::Real=1.15,
        srp_Am::Real=0.002,
        srp_P::Real=4.56e-6;
        naif_frame::String="J2000",
        abcorr::String="NONE",
        AU = 149597870.7,
        use_canonical::Bool = true,
    )
        if use_canonical
            # scaled mus
            mus_scaled = mus / mus[1]

            # calculate tstar and vstar
            vstar = sqrt(mus[1]/lstar)
            tstar = lstar/vstar
        else
            mus_scaled = mus
            lstar = 1.0
            vstar = 1.0
            tstar = 1.0
        end

        # compute magnitude scalar for SRP
        k_srp = (AU/lstar)^2 * (srp_P * srp_cr * srp_Am / 1000) * (tstar^2/lstar)

        # initialize object
        new(
            et0,
            lstar,
            tstar,
            vstar,
            mus,
            mus_scaled,
            k_srp,
            naif_ids,
            naif_frame,
            abcorr,
            nothing,
            zeros(MVector{3 * (length(mus)-1),Float64})
        )
    end
end
