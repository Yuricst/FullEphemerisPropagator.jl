"""
Equations of motion
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

    # constructor
    function Nbody_params(
        et0,
        lstar::Real,
        mus::Vector{Float64},
        naif_ids::Vector{String};
        naif_frame::String="J2000",
        abcorr::String="NONE",
    )
        # scaled mus
        mus_scaled = mus / mus[1]
        # calculate tstar and vstar
        vstar = sqrt(lstar/mus[1])
        tstar = lstar/vstar
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
        )
    end
end


"""
Compute third-body acceleration via Battin's formula
"""
function third_body_accel(r_spacecraft, r_3body, mu_3body)
    s = r_spacecraft - r_3body
    q = dot(r_spacecraft, r_spacecraft - 2s)/dot(s, s)
    F = q * (3 + 3q + q^2)/(1 + sqrt(1+q)^3)
    return -mu_3body/norm(r_3body)^3 * (r_spacecraft + F*s)
end


"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
"""
function eom_Nbody_SPICE!(du, u, params, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # position derivatives
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # velocity derivatives
    du[4] = -mu_r3 * u[1]
    du[5] = -mu_r3 * u[2]
    du[6] = -mu_r3 * u[3]

    # third-body effects
    for i = 2:length(params.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            params.naif_ids[i],
            params.et0 + t*params.tstar,
            params.naif_frame,
            params.abcorr,
            params.naif_ids[1]
        )
        pos_3body /= params.lstar   # re-scale

        # compute third-body perturbation
        du[4:6] += third_body_accel(u[1:3], pos_3body, params.mus_scaled[i])
    end
end


"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
This function propagates the concatenated state and STM.
"""
function eom_NbodySTM_SPICE!(du, u, params, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # position derivatives
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # velocity derivatives
    du[4] = -mu_r3 * u[1]
    du[5] = -mu_r3 * u[2]
    du[6] = -mu_r3 * u[3]

    # third-body effects
    Rs = zeros(3, length(params.mus_scaled)-1)
    for i = 2:length(params.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            params.naif_ids[i],
            params.et0 + t*params.tstar,
            params.naif_frame,
            params.abcorr,
            params.naif_ids[1]
        )
        pos_3body /= params.lstar   # re-scale
        Rs[:,i-1] = pos_3body
        
        # compute third-body perturbation
        du[4:6] += third_body_accel(u[1:3], pos_3body, params.mus_scaled[i])
    end

    # stm derivatives
    jacobian = params.f_jacobian([u[1:3]..., params.mus_scaled..., Rs])
    du[7:42] = reshape(jacobian * reshape(u[7:42], (6,6)), 36)
end
