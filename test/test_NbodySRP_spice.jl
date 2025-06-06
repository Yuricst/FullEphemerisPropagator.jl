"""
Test integrating N-body dynamics with SPICE call within eom
"""

if !@isdefined(FullEphemerisPropagator)
    include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))
end

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
        use_srp = true,
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
        -0.9707266474653355
        -0.5034448179024911
        -1.775070433549609
         0.1292918188529246
        -0.5984316125180164
         0.25702901382064053
    ]
    @test norm(sol.u[end] - xf_check) < 1e-11
end

test_Nbody_spice()
