"""
Functions to generate symbolic Jacobian
"""


function symbolic_Nbody_jacobian(N::Int)
    # define symbolic variables
    Symbolics.@variables x y z #state[1:3]
    Symbolics.@variables mus[1:N]
    Symbolics.@variables Rs[1:3, 1:N-1]

    # define accelerations
    rvec = [x,y,z]
    rnorm = norm(rvec)
    ax = -mus[1]*x/rnorm^3
    ay = -mus[1]*y/rnorm^3
    az = -mus[1]*z/rnorm^3
    
    # third-body accelerations
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


function symbolic_NbodySRP_jacobian(N::Int)
    # define symbolic variables
    Symbolics.@variables x y z #state[1:3]
    Symbolics.@variables mus[1:N]
    Symbolics.@variables Rs[1:3, 1:N-1]
    Symbolics.@variables R_sun[1:3]
    Symbolics.@variables k_srp

    # define accelerations
    rvec = [x,y,z]
    rnorm = norm(rvec)
    ax = -mus[1]*x/rnorm^3
    ay = -mus[1]*y/rnorm^3
    az = -mus[1]*z/rnorm^3
    
    # third-body accelerations
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

    # SRP acceleration
    x_Sun2sc = x - R_sun[1]
    y_Sun2sc = y - R_sun[2]
    z_Sun2sc = z - R_sun[3]
    r_relative_norm3 = (x_Sun2sc^2 + y_Sun2sc^2 + z_Sun2sc^2)^(3/2)
    ax += k_srp * x_Sun2sc / r_relative_norm3
    ay += k_srp * y_Sun2sc / r_relative_norm3
    az += k_srp * z_Sun2sc / r_relative_norm3
    
    jac = [
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        Symbolics.derivative(ax, x) Symbolics.derivative(ax, y) Symbolics.derivative(ax, z) 0 0 0;
        Symbolics.derivative(ay, x) Symbolics.derivative(ay, y) Symbolics.derivative(ay, z) 0 0 0;
        Symbolics.derivative(az, x) Symbolics.derivative(az, y) Symbolics.derivative(az, z) 0 0 0;
    ]

    arguments = [x, y, z, mus..., Rs, R_sun, k_srp]
    f_jacobian, _ = Symbolics.build_function(jac, arguments; expression = Val{false})
    return Symbolics.eval(f_jacobian)
end