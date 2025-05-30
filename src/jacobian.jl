"""Define Jacobioan fetching functions"""


"""
Compute Jacobian from dynamics in propagator model.
Evaluated at epoch `et0 + t*tstar` at state `u`
"""
function jacobian(
    propagator::PropagatorSTM,
    et0::Float64,
    t::Real,
    u::Vector,
)
    # get third-body positions
    for i = 2:length(propagator.parameters.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            propagator.parameters.naif_ids[i],
            et0 + t * propagator.parameters.tstar,
            propagator.parameters.naif_frame,
            propagator.parameters.abcorr,
            propagator.parameters.naif_ids[1]
        )
        pos_3body /= propagator.parameters.lstar                       # re-scale
        propagator.parameters.Rs[1+3(i-2):3(i-1)] .= pos_3body
    end

    Uxx = propagator.parameters.f_jacobian(
        u[1:3]...,
        propagator.parameters.mus_scaled...,
        propagator.parameters.Rs...)
    return vcat(hcat(zeros(3,3), I(3)), hcat(reshape(Uxx, (3,3))', zeros(3,3)))
end


function jacobian(
    parameters::Nbody_params,
    et0::Float64,
    t::Real,
    u::Vector,
)
    # get third-body positions
    for i = 2:length(parameters.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            parameters.naif_ids[i],
            et0 + t * parameters.tstar,
            parameters.naif_frame,
            parameters.abcorr,
            parameters.naif_ids[1]
        )
        pos_3body /= parameters.lstar                       # re-scale
        parameters.Rs[1+3(i-2):3(i-1)] .= pos_3body
    end

    Uxx = parameters.f_jacobian(
        u[1:3]...,
        parameters.mus_scaled...,
        parameters.Rs...)
    return vcat(hcat(zeros(3,3), I(3)), hcat(reshape(Uxx, (3,3))', zeros(3,3)))
end


function jacobian(
    parameters::NbodySRP_params,
    et0::Float64,
    t::Real,
    u::Vector,
)
    # get third-body positions
    R_sun = zeros(3)
    for i = 2:length(parameters.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            parameters.naif_ids[i],
            et0 + t * parameters.tstar,
            parameters.naif_frame,
            parameters.abcorr,
            parameters.naif_ids[1]
        )
        pos_3body /= parameters.lstar                       # re-scale
        parameters.Rs[1+3(i-2):3(i-1)] .= pos_3body
        if parameters.naif_ids[i] == "10"
            R_sun[:] = pos_3body
        end
    end

    Uxx = parameters.f_jacobian(
        u[1:3]...,
        parameters.mus_scaled...,
        parameters.Rs...,
        R_sun...,
        parameters.k_srp)
    return vcat(hcat(zeros(3,3), I(3)), hcat(reshape(Uxx, (3,3))', zeros(3,3)))
end