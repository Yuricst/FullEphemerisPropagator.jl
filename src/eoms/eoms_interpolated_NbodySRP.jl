"""Equations of motion for N-body problem with solar radiation pressure (SRP)"""


"""
N-body equations of motion with SRP, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
"""
function eom_NbodySRP_SPICE(u, params::InterpolatedNbodySRP_params, t)
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

        # add SRP
        if params.naif_ids[i] == "10"
            r_relative = u[1:3] - r_3body   # Sun -> spacecraft vector
            du[4:6] += params.k_srp * r_relative/norm(r_relative)^3
            dvx += params.k_srp * r_relative[1]/norm(r_relative)^3
            dvy += params.k_srp * r_relative[2]/norm(r_relative)^3
            dvz += params.k_srp * r_relative[3]/norm(r_relative)^3
        end
    end
    return SA[dx,dy,dz,dvx,dvy,dvz]
end



"""
N-body equations of motion with SRP, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
"""
function eom_NbodySRP_SPICE!(du, u, params::InterpolatedNbodySRP_params, t)
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

        # add SRP
        if params.naif_ids[i] == "10"
            r_relative = u[1:3] - r_3body   # Sun -> spacecraft vector
            du[4:6] += params.k_srp * r_relative/norm(r_relative)^3
        end
    end
    return nothing
end


"""
N-body equations of motion with SRP, using SPICE query for third-body positions.
This function signature is compatible with `DifferentialEquations.jl`.
This function propagates the concatenated state and STM.
"""
function eom_NbodySRP_STM_SPICE!(du, u, params::InterpolatedNbodySRP_params, t)
    # compute coefficient
    mu_r3 = (params.mus_scaled[1] / norm(u[1:3])^3)

    # position derivatives
    du[1:3] .= u[4:6]

    # velocity derivatives
    du[4:6] .= -mu_r3 * u[1:3]

    # third-body effects
    R_sun = zeros(3)
    for i = 2:length(params.mus_scaled)
        # get position of third body
        r_3body = get_pos(params.ephem_dict[i], params.et0 + t*params.tstar)
        params.Rs[1+3(i-2):3(i-1)] .= r_3body
        
        # compute third-body perturbation
        du[4:6] += third_body_accel(u[1:3], r_3body, params.mus_scaled[i])

        # add SRP
        if params.naif_ids[i] == "10"
            r_relative = u[1:3] - r_3body   # Sun -> spacecraft vector
            du[4:6] += params.k_srp * r_relative/norm(r_relative)^3
            R_sun .= r_3body
        end
    end

    # stm derivatives
    Uxx = params.f_jacobian(u[1:3]..., params.mus_scaled..., params.Rs..., R_sun..., params.k_srp)
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