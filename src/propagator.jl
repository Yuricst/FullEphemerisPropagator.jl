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
        naif_ids::Vector{String},
        srp_cr::Real=1.15,
        srp_Am::Real=0.002,
        srp_P::Real=4.56e-6;
        use_srp::Bool=false,
        naif_frame::String="J2000",
        abcorr::String="NONE",
        reltol::Float64=1e-12,
        abstol::Float64=1e-12,
    )
        @assert length(mus) == length(naif_ids) "Length of mus and naif_ids should be equal!"

        # construct parameters
        if use_srp == false
            eom! = eom_Nbody_SPICE!
            parameters = Nbody_params(
                str2et("2000-01-01T00:00:00"),
                lstar,
                mus,
                naif_ids;
                naif_frame=naif_frame,
                abcorr=abcorr
            )
        else
            eom! = eom_NbodySRP_SPICE!
            parameters = NbodySRP_params(
                str2et("2000-01-01T00:00:00"),
                lstar,
                mus,
                naif_ids,
                srp_cr,
                srp_Am,
                srp_P;
                naif_frame=naif_frame,
                abcorr=abcorr
            )
        end

        # construct ODE problem
        problem = ODEProblem(
            eom!,
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
Propagator struct holds ODEProblem along with solve parameters.
"""
mutable struct PropagatorSTM <: FullEphemPropagator
    problem
    parameters::FullEphemParameters
    method
    reltol::Float64
    abstol::Float64
    
    # constructor
    function PropagatorSTM(
        method,
        lstar::Real,
        mus::Vector{Float64},
        naif_ids::Vector{String};
        naif_frame::String="J2000",
        abcorr::String="NONE",
        reltol::Float64=1e-12,
        abstol::Float64=1e-12,
    )
        @assert length(mus) == length(naif_ids) "Length of mus and naif_ids should be equal!"

        # construct parameters
        parameters = Nbody_params(
            str2et("2000-01-01T00:00:00"),
            lstar,
            mus,
            naif_ids;
            naif_frame=naif_frame,
            abcorr=abcorr
        )

        # mutate parameter `f_jacobian`
        parameters.f_jacobian = symbolic_Nbody_jacobian(length(mus))

        # construct ODE problem
        problem = ODEProblem(
            eom_NbodySTM_SPICE!,
            vcat([1.0, 0.0, 0.0, 0.5, 1.0, 0.0], ones(36)),  # placeholder for u0
            [0.0, 1.0],                                      # placeholder for tspan
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


function pretty(propagator::Union{Propagator,PropagatorSTM})
    println("Full-ephemeris integrator")
    println("   naif_ids     : ", propagator.parameters.naif_ids)
    @printf("   canonical LU : %.1e\n", propagator.parameters.lstar)
    @printf("   canonical TU : %.1e\n", propagator.parameters.tstar)
    println("   method       : ", propagator.method)
    @printf("   reltol       : %.1e\n", propagator.reltol)
    @printf("   abstol       : %.1e\n", propagator.abstol)
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.
"""
function propagate(
    propagator::Propagator,
    et0::Float64,
    tspan::Tuple{Real,Real},
    u0::Vector;
    callback = nothing,
    kwargs...
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    # modify value for initial epoch
    propagator.parameters.et0 = et0

    # mutate ODEProblem 
    modified_problem = remake(
        propagator.problem;
        u0 = u0,
        tspan = tspan,
        p = propagator.parameters
    )

    # solve and return results
    return solve(modified_problem,
                 propagator.method;
                 callback = callback,
                 reltol = propagator.reltol,
                 abstol = propagator.abstol,
                 kwargs...)
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.
"""
function propagate(
    propagator::PropagatorSTM,
    et0::Float64,
    tspan::Tuple{Real,Real},
    u0::Vector;
    callback = nothing,
    stm0 = nothing,
    kwargs...
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    if isnothing(stm0) == true
        u0_stm = vcat(u0, reshape(I(6), 36))
    else
        u0_stm = vcat(u0, reshape(stm0, 36))
    end

    # modify value for initial epoch
    propagator.parameters.et0 = et0

    # mutate ODEProblem 
    modified_problem = remake(
        propagator.problem;
        u0 = u0_stm,
        tspan = tspan,
        p = propagator.parameters
    )

    # solve and return results
    return solve(modified_problem,
                 propagator.method;
                 callback = callback,
                 reltol = propagator.reltol,
                 abstol = propagator.abstol,
                 kwargs...)
end
