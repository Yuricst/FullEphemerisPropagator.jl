"""
Test integrating N-body dynamics with SPICE call within eom.
Uses low-level API
"""

test_Nbody_spice_lowapi = function()
    # define parameters
    mus = [
        4.9028000661637961E+03,
        3.9860043543609598E+05,
        1.3271244004193938E+11,
    ]
    naif_ids = ["301", "399", "10"]
    naif_frame = "J2000"
    abcorr = "NONE"
    lstar = 3000.0

    et0 = str2et("2020-01-01T00:00:00")
    parameters = FullEphemerisPropagator.Nbody_params(
        et0,
        lstar,
        mus,
        naif_ids;
        naif_frame=naif_frame,
        abcorr=abcorr
    )

    # initial state (in canonical scale)
    u0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]

    # time span (in canonical scale)
    tspan = (0.0, 7*86400/parameters.tstar)

    # solve
    prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, parameters)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    u_check = [
        0.522315099783345
        2.096145050620713
       -0.16366006795323232
       -0.40936130082155736
        0.25386267312152233
       -0.16005652013568536
    ]
    @test norm(sol.u[end] - u_check) < 1e-12
end

test_Nbody_spice_lowapi()