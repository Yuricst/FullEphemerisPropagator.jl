"""
Test integrating N-body dynamics with SPICE call within eom
"""

test_Nbody_spice = function()

    # define parameters
    naif_ids = ["301", "399", "10"]
    mus = [bodvrd(ID, "GM", 1)[1] for ID in naif_ids]
    naif_frame = "J2000"
    abcorr = "NONE"
    lstar = 3000.0

    # instantiate propagator
    prop = FullEphemerisPropagator.Propagator(
        Vern7(),
        lstar,
        mus,
        naif_ids;
        naif_frame = naif_frame,
        reltol = 1e-12,
        abstol = 1e-12,
    )

    # initial epoch
    et0 = str2et("2020-01-01T00:00:00")

    # initial state (in canonical scale)
    u0_dim = [2200.0, 0.0, 4200.0, 0.03, 1.1, 0.1]
    u0 = FullEphemerisPropagator.dim2nondim(prop, u0_dim)

    # time span (in canonical scale)
    tspan = (0.0, 30*86400/prop.parameters.tstar)

    # solve
    tevals = LinRange(tspan[1], tspan[2], 15000)   # optionally specify when to query states
    sol = FullEphemerisPropagator.propagate(prop, et0, tspan, u0; saveat=tevals)
    xf_check = [
        -0.9706832470466424
        -0.5034938627169759
        -1.775010908894415
         0.12930004488210714
        -0.5984476217012531
         0.25704051772397823
    ]
    @test norm(sol.u[end] - xf_check) < 1e-11
end

test_Nbody_spice()
