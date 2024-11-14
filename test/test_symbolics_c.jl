using Symbolics


# function lotka_volterra!(du, u, p, t)
#   x, y = u
#   α, β, δ, γ = p
#   du[1] = dx = α*x - β*x*y
#   du[2] = dy = -δ*y + γ*x*y
# end

# @variables t du[1:2] u[1:2] p[1:4]
# du = collect(du)
# lotka_volterra!(du, u, p, t)
# du


# f = build_function(du, u, p, t, target=Symbolics.CTarget(), expression=Val{false})

# # test
# du = rand(2); du2 = rand(2)
# u = rand(2)
# p = rand(4)
# t = rand()
# f(du, u, p, t)
# lotka_volterra!(du2, u, p, t)
# du == du2 # true!


function get_func(N::Int)
    # define symbolic variables
    Symbolics.@variables x y z #state[1:3]
    Symbolics.@variables mus[1:N]
    Symbolics.@variables Rs[3*(N-1)]  #[1:3, 1:N-1]
    
    # define accelerations
    rnorm = sqrt(x^2 + y^2 + z^2)
    ax = -mus[1]*x/rnorm^3
    ay = -mus[1]*y/rnorm^3
    az = -mus[1]*z/rnorm^3

    # third-body accelerations
    for i = 2:N
        R3 = sqrt(Rs[1+3(i-2)]^2 + Rs[2+3(i-2)]^2 + Rs[3+3(i-2)]^2)^3
        # s = rvec - Rs[:,i-1]
        snorm2 = (x - Rs[1+3(i-2)])^2 + (y - Rs[2+3(i-2)])^2 + (z - Rs[3+3(i-2)])^2

        #rvec_2s = [x,y,z] - 2*s
        #q = (x*rvec_2s[1] + y*rvec_2s[2] + z*rvec_2s[3]) / snorm2

        q = (x * (x - 2*(x - Rs[1+3(i-2)])) +
            y * (y - 2*(y - Rs[2+3(i-2)])) +
            z * (z - 2*(z - Rs[3+3(i-2)]))) / snorm2
        F = q * (3 + 3q + q^2)/(1 + sqrt(1+q)^3)

        ax += -mus[i] / R3 * (x + F*(x - Rs[1+3(i-2)]))
        ay += -mus[i] / R3 * (y + F*(y - Rs[2+3(i-2)]))
        az += -mus[i] / R3 * (z + F*(z - Rs[3+3(i-2)]))
        ax += snorm2
    end

    Uxx = [
        Symbolics.derivative(ax, x) Symbolics.derivative(ax, y) Symbolics.derivative(ax, z);
        Symbolics.derivative(ay, x) Symbolics.derivative(ay, y) Symbolics.derivative(ay, z);
        Symbolics.derivative(az, x) Symbolics.derivative(az, y) Symbolics.derivative(az, z);
    ]

    arguments = [x, y, z, mus..., Rs...]
    return Symbolics.build_function(Uxx, arguments;
                                    target=Symbolics.CTarget(),
                                    expression=Val{false})
end

get_func(3)