"""
Integrator object
"""

abstract type FullEphemPropagator end


"""
Propagator struct holds ODEProblem along with solve parameters.
"""
mutable struct Propagator <: FullEphemPropagator
    problem
    parameters::FullEphemParameters
    method
    reltol::Float64
    abstol::Float64
    
    # constructor
    function Propagator(
        method,
        lstar::Real,
        mus::Vector{Float64},
        naif_ids::Vector{String};
        naif_frame::String="J2000",
        abcorr::String="NONE",
        reltol::Float64=1e-12,
        abstol::Float64=1e-12,
    )
        # construct parameters
        parameters = FullEphemerisPropagator.Nbody_params(
            str2et("2000-01-01T00:00:00"),
            lstar,
            mus,
            naif_ids;
            naif_frame=naif_frame,
            abcorr=abcorr
        )

        # construct ODE problem
        problem = ODEProblem(
            eom_Nbody_SPICE!,
            [1.0, 0.0, 0.0, 0.5, 1.0, 0.0],  # placeholder for u0
            [0.0, 1.0],                      # placeholder for tspan
            parameters,
        )

        # initialize
        new(problem,
            parameters,
            method,
            reltol,
            abstol,
        )
    end
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.
"""
function propagate(
    propagator::FullEphemPropagator,
    et0::Float64,
    u0::Vector,
    tspan::Tuple{Real,Real};
    callback = nothing,
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    # modify value for initial epoch
    propagator.parameters.et0 = et0

    # mutate ODEProblem 
    modified_problem = remake(propagator.problem;
        u0 = u0,
        tspan = tspan,
        p = propagator.parameters
    )

    # solve and return results
    return solve(modified_problem, propagator.method;
                 callback = callback,
                 reltol=propagator.reltol, abstol=propagator.abstol)
end