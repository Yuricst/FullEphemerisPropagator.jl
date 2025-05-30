"""Equations of motion for N-body problem with interpolation of third-body positions"""


"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
This is a static version of the function.
"""
function eom_Nbody_SPICE(u, params::InterpolatedNbodyParams, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # position derivatives
    dx, dy, dz = u[4:6]

    # velocity derivatives
    dvx, dvy, dvz = -mu_r3 * u[1:3]

    # third-body effects
    for i = 2:length(params.mus_scaled)
        # get position of third body
        r_3body = get_pos(params.ephem_dict[i], params.et0 + t*params.tstar)
        
        # compute third-body perturbation
        a_3bd = third_body_accel(u[1:3], r_3body, params.mus_scaled[i])
        dvx += a_3bd[1]
        dvy += a_3bd[2]
        dvz += a_3bd[3]
    end
    return SA[dx,dy,dz,dvx,dvy,dvz]
end



"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
"""
function eom_Nbody_SPICE!(du, u, params::InterpolatedNbodyParams, t)
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
        r_3body = get_pos(params.ephem_dict[i], params.et0 + t*params.tstar)

        # compute third-body perturbation
        du[4:6] += third_body_accel(u[1:3], r_3body, params.mus_scaled[i])
    end

    return nothing
end



"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
This function propagates the concatenated state and STM.
"""
function eom_Nbody_STM_SPICE!(du, u, params::InterpolatedNbodyParams, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # position derivatives
    du[1:3] .= u[4:6]

    # velocity derivatives
    du[4:6] .= -mu_r3 * u[1:3]

    # third-body effects
    #Rs = zeros(3, length(params.mus_scaled)-1)
    for i = 2:length(params.mus_scaled)
        # get position of third body
        r_3body = get_pos(params.ephem_dict[i], params.et0 + t*params.tstar)
        params.Rs[1+3(i-2):3(i-1)] .= r_3body
        
        # compute third-body perturbation
        du[4:6] += third_body_accel(u[1:3], r_3body, params.mus_scaled[i])
    end

    # stm derivatives
    Uxx = params.f_jacobian(u[1:3]..., params.mus_scaled..., params.Rs...)
    du[7:12]  .= u[25:30]
    du[13:18] .= u[31:36]
    du[19:24] .= u[37:42]
    
    du[25] = Uxx[1]*u[7]  + Uxx[2]*u[13] + Uxx[3]*u[19]
    du[26] = Uxx[1]*u[8]  + Uxx[2]*u[14] + Uxx[3]*u[20]
    du[27] = Uxx[1]*u[9]  + Uxx[2]*u[15] + Uxx[3]*u[21]
    du[28] = Uxx[1]*u[10] + Uxx[2]*u[16] + Uxx[3]*u[22]
    du[29] = Uxx[1]*u[11] + Uxx[2]*u[17] + Uxx[3]*u[23]
    du[30] = Uxx[1]*u[12] + Uxx[2]*u[18] + Uxx[3]*u[24]
    
    du[31] = Uxx[4]*u[7]  + Uxx[5]*u[13] + Uxx[6]*u[19]
    du[32] = Uxx[4]*u[8]  + Uxx[5]*u[14] + Uxx[6]*u[20]
    du[33] = Uxx[4]*u[9]  + Uxx[5]*u[15] + Uxx[6]*u[21]
    du[34] = Uxx[4]*u[10] + Uxx[5]*u[16] + Uxx[6]*u[22]
    du[35] = Uxx[4]*u[11] + Uxx[5]*u[17] + Uxx[6]*u[23]
    du[36] = Uxx[4]*u[12] + Uxx[5]*u[18] + Uxx[6]*u[24]

    du[37] = Uxx[7]*u[7]  + Uxx[8]*u[13] + Uxx[9]*u[19]
    du[38] = Uxx[7]*u[8]  + Uxx[8]*u[14] + Uxx[9]*u[20]
    du[39] = Uxx[7]*u[9]  + Uxx[8]*u[15] + Uxx[9]*u[21]
    du[40] = Uxx[7]*u[10] + Uxx[8]*u[16] + Uxx[9]*u[22]
    du[41] = Uxx[7]*u[11] + Uxx[8]*u[17] + Uxx[9]*u[23]
    du[42] = Uxx[7]*u[12] + Uxx[8]*u[18] + Uxx[9]*u[24]
    return nothing
end



"""
N-body equations of motion, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
This function propagates the concatenated state and STM.
This is a static version of the function.
"""
function eom_Nbody_STM_SPICE(u, params::InterpolatedNbodyParams, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # velocity derivatives
    dvx, dvy, dvz = -mu_r3 * u[1:3]

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
        params.Rs[1+3(i-2):3(i-1)] .= pos_3body
        
        # compute third-body perturbation
        a_3bd = third_body_accel(u[1:3], pos_3body, params.mus_scaled[i])
        dvx += a_3bd[1]
        dvy += a_3bd[2]
        dvz += a_3bd[3]
    end

    # stm derivatives
    Uxx = params.f_jacobian(u[1:3]..., params.mus_scaled..., params.Rs...)

    return SA[u[4:6]..., dvx, dvy, dvz,
              u[25:30]..., u[31:36]..., u[37:42]...,
              Uxx[1]*u[7]  + Uxx[2]*u[13] + Uxx[3]*u[19],
              Uxx[1]*u[8]  + Uxx[2]*u[14] + Uxx[3]*u[20],
              Uxx[1]*u[9]  + Uxx[2]*u[15] + Uxx[3]*u[21],
              Uxx[1]*u[10] + Uxx[2]*u[16] + Uxx[3]*u[22],
              Uxx[1]*u[11] + Uxx[2]*u[17] + Uxx[3]*u[23],
              Uxx[1]*u[12] + Uxx[2]*u[18] + Uxx[3]*u[24],
              Uxx[4]*u[7]  + Uxx[5]*u[13] + Uxx[6]*u[19],
              Uxx[4]*u[8]  + Uxx[5]*u[14] + Uxx[6]*u[20],
              Uxx[4]*u[9]  + Uxx[5]*u[15] + Uxx[6]*u[21],
              Uxx[4]*u[10] + Uxx[5]*u[16] + Uxx[6]*u[22],
              Uxx[4]*u[11] + Uxx[5]*u[17] + Uxx[6]*u[23],
              Uxx[4]*u[12] + Uxx[5]*u[18] + Uxx[6]*u[24],
              Uxx[7]*u[7]  + Uxx[8]*u[13] + Uxx[9]*u[19],
              Uxx[7]*u[8]  + Uxx[8]*u[14] + Uxx[9]*u[20],
              Uxx[7]*u[9]  + Uxx[8]*u[15] + Uxx[9]*u[21],
              Uxx[7]*u[10] + Uxx[8]*u[16] + Uxx[9]*u[22],
              Uxx[7]*u[11] + Uxx[8]*u[17] + Uxx[9]*u[23],
              Uxx[7]*u[12] + Uxx[8]*u[18] + Uxx[9]*u[24],
            ]
end