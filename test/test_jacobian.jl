"""
Test jacobian function
"""

if !@isdefined(FullEphemerisPropagator)
    include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))
end

test_jacobian = function ()

    naif_frame = "J2000"
    lstar = 3000
    mus = [
        4.9028000661637961E+03,
        3.9860043543609598E+05,
        1.3271244004193938E+11,
    ]
    naif_ids = ["301", "399", "10"]
    N = length(naif_ids)

    # integrate with STM
    prop = FullEphemerisPropagator.PropagatorSTM(
        Vern9(),
        lstar,
        mus,
        naif_ids;
        naif_frame = naif_frame,
        reltol = 1e-12,
        abstol = 1e-12,
    )

    # initial epoch
    et0 = 840463475.5117936

    # initial state (in canonical scale)
    # u0_dim = [3200.0, 0.0, 4200.0, 0.0, 0.77, 0.0]
    # u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)
    u0 = [
        -2.5019204591096096,
        14.709398066624694,
        -18.59744250295792,
        5.62688812721852e-2,
        1.439926311669468e-2,
        3.808273517470642e-3
    ]

    # time span (in canonical scale)
    tspan = (0.0, 6.55 * 5 *86400/prop.parameters.tstar)

    # evaluate jacobian at initial time
    jac = FullEphemerisPropagator.jacobian(prop, et0, 0.0, u0)

    Uxx_check = [
        -7.22461e-5   2.03104e-5   4.94217e-5;
        2.03104e-5   1.45477e-5  -7.25469e-5;
        4.94217e-5  -7.25469e-5   5.76984e-5;
    ]
    
    @test jac[1:3,1:3] == zeros(3,3)
    @test jac[1:3,4:6] == I(3)
    @test norm(jac[4:6,1:3] - Uxx_check, Inf) < 1e-10
    @test jac[4:6,4:6] == zeros(3,3)
end

test_jacobian()