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

        # # compute radius to third body
        du[4:6] += third_body_accel(u[1:3], pos_3body, params.mus_scaled[i])
    end
end


function symbolic_Nbody_jacobian(N::Int)
    # define symbolic variables
    Symbolics.@variables x y z #state[1:3]
    Symbolics.@variables mus[1:N]
    Symbolics.@variables Rs[1:3, 1:N-1]

    # unpack state
    #x, y, z = state

    # define accelerations
    rvec = [x,y,z]
    rnorm = norm(rvec)
    ax = -mus[1]*x/rnorm^3
    ay = -mus[1]*y/rnorm^3
    az = -mus[1]*z/rnorm^3
    
    for i = 2:N
        R3 = sqrt(Rs[1,i-1]^2 + Rs[2,i-1]^2 + Rs[3,i-1]^2)^3
        # s = rvec - Rs[:,i-1]
        snorm2 = (x - Rs[1,i-1])^2 + (y - Rs[2,i-1])^2 + (z - Rs[3,i-1])^2

        #rvec_2s = [x,y,z] - 2*s
        #q = (x*rvec_2s[1] + y*rvec_2s[2] + z*rvec_2s[3]) / snorm2

        q = (x * (x - 2*(x - Rs[1,i-1])) +
             y * (y - 2*(y - Rs[2,i-1])) +
             z * (z - 2*(z - Rs[3,i-1]))) / snorm2
        F = q * (3 + 3q + q^2)/(1 + sqrt(1+q)^3)

        ax += -mus[i] / R3 * (x + F*(x - Rs[1,i-1]))
        ay += -mus[i] / R3 * (y + F*(y - Rs[2,i-1]))
        az += -mus[i] / R3 * (z + F*(z - Rs[3,i-1]))
    end
    
    jac = [
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        Symbolics.derivative(ax, x) Symbolics.derivative(ax, y) Symbolics.derivative(ax, z) 0 0 0;
        Symbolics.derivative(ay, x) Symbolics.derivative(ay, y) Symbolics.derivative(ay, z) 0 0 0;
        Symbolics.derivative(az, x) Symbolics.derivative(az, y) Symbolics.derivative(az, z) 0 0 0;
    ]

    arguments = [x, y, z, mus..., Rs]
    f_jacobian, _ = Symbolics.build_function(jac, arguments; expression = Val{false})
    return Symbolics.eval(f_jacobian)
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
        
        # # compute radius to third body
        du[4:6] += third_body_accel(u[1:3], pos_3body, params.mus_scaled[i])
    end

    # stm derivatives
    jacobian = params.f_jacobian([u[1:3]..., params.mus_scaled..., Rs])
    du[7:42] = reshape(jacobian * reshape(u[7:42], (6,6)), 36)
end
