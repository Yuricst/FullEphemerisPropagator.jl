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
    use_srp::Bool
    use_sa::Bool

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
        use_sa::Bool=false,
    )
        @assert length(mus) == length(naif_ids) "Length of mus and naif_ids should be equal!"

        # construct parameters
        if use_srp == false
            if use_sa == true
                eom! = eom_Nbody_SPICE
            else
                eom! = eom_Nbody_SPICE!
            end
            parameters = Nbody_params(
                str2et("2000-01-01T00:00:00"),
                lstar,
                mus,
                naif_ids;
                naif_frame=naif_frame,
                abcorr=abcorr
            )
        else
            if use_sa == true
                eom! = eom_NbodySRP_SPICE
            else
                eom! = eom_NbodySRP_SPICE!
            end
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
            use_srp,
            use_sa,
        )
    end
end


"""
Overload method for showing Propagator
"""
function Base.show(io::IO, propagator::Propagator)
    println("Propagator struct")
    @printf("    method  : %s\n", string(propagator.method))
    @printf("    reltol  : %1.4f\n", propagators.reltol)
    @printf("    abstol  : %1.4f\n", propagators.abstol)
    @printf("    use_srp : %s\n", string(propagator.use_srp))
    @printf("    use_sa  : %s\n", string(propagator.use_sa))
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
    use_srp::Bool
    use_sa::Bool
    
    # constructor
    function PropagatorSTM(
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
        use_sa::Bool=false,
    )
        @assert length(mus) == length(naif_ids) "Length of mus and naif_ids should be equal!"

        # construct parameters
        if use_srp == false
            if use_sa == true
                eom! = eom_Nbody_STM_SPICE
            else
                eom! = eom_Nbody_STM_SPICE!
            end
            parameters = Nbody_params(
                str2et("2000-01-01T00:00:00"),
                lstar,
                mus,
                naif_ids;
                naif_frame=naif_frame,
                abcorr=abcorr
            )
            parameters.f_jacobian = symbolic_Nbody_jacobian(length(mus))
        else
            if use_sa == true
                error("Not implemented yet!")
                #eom! = eom_Nbody_SPICE
            else
                eom! = eom_NbodySRP_STM_SPICE!
            end
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
            parameters.f_jacobian = symbolic_NbodySRP_jacobian(length(mus))
        end

        # construct ODE problem
        problem = ODEProblem(
            eom!,
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
            use_srp,
            use_sa,
        )
    end
end


"""
Overload method for showing PropagatorSTM
"""
function Base.show(io::IO, propagator::PropagatorSTM)
    println("PropagatorSTM struct")
    @printf("    method  : %s\n", string(propagator.method))
    @printf("    reltol  : %1.4f\n", propagators.reltol)
    @printf("    abstol  : %1.4f\n", propagators.abstol)
    @printf("    use_srp : %s\n", string(propagator.use_srp))
    @printf("    use_sa  : %s\n", string(propagator.use_sa))
end


"""
Propagate initial state `u0` from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.

Additional keyworded arguments for DifferentialEquations.solve() can be passed.
See: https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple%7BSciMLBase.AbstractDEProblem,%20Vararg%7BAny%7D%7D

# Arguments
- `propagator::Propagator`: Propagator object
- `et0::Float64`: Initial epoch in ephemeris time, in seconds
- `tspan::Tuple{Real,Real}`: Time span to propagate, in canonical time units
- `u0::Vector`: Initial state vector
- `callback::Union{Nothing,Function}`: Optional callback function
- `kwargs...`: Additional keyworded arguments for DifferentialEquations.solve()
"""
function propagate(
    propagator::Propagator,
    et0::Float64,
    tspan::Tuple{Real,Real},
    u0::Union{Vector,SVector{6,Float64}};
    callback = nothing,
    kwargs...
)
    @assert length(u0) == 6 "u0 should be length 6, not $(length(u0))!"
    @assert length(tspan) == 2 "tspan should be length 2, not $(length(tspan))!"

    # modify value for initial epoch
    propagator.parameters.et0 = et0

    # mutate ODEProblem 
    #u0 = @MVector [el for el in u0]
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
Propagate initial state `u0` and STM from `tspan[1]` to `tspan[2]`.
The initial state should be given as `u0 = [x,y,z,vx,vy,vz]`.

Additional keyworded arguments for DifferentialEquations.solve() can be passed.
See: https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple%7BSciMLBase.AbstractDEProblem,%20Vararg%7BAny%7D%7D    

# Arguments
- `propagator::Propagator`: Propagator object
- `et0::Float64`: Initial epoch in ephemeris time, in seconds
- `tspan::Tuple{Real,Real}`: Time span to propagate, in canonical time units
- `u0::Vector`: Initial state vector
- `callback::Union{Nothing,Function}`: Optional callback function
- `kwargs...`: Additional keyworded arguments for DifferentialEquations.solve()
"""
function propagate(
    propagator::PropagatorSTM,
    et0::Float64,
    tspan::Tuple{Real,Real},
    u0::Union{Vector,SVector{6,Float64}};
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

    if propagator.use_sa
        u0_stm = SA[u0_stm...]
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



"""
Compute Jacobian from dynamics in propagator model.
Evaluated at epoch `et0 + t*tstar` at state `u`
"""
function jacobian(
    propagator::PropagatorSTM,
    et0::Float64,
    t::Real,
    u::Vector,
)
    # get third-body positions
    for i = 2:length(propagator.parameters.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            propagator.parameters.naif_ids[i],
            et0 + t * propagator.parameters.tstar,
            propagator.parameters.naif_frame,
            propagator.parameters.abcorr,
            propagator.parameters.naif_ids[1]
        )
        pos_3body /= propagator.parameters.lstar                       # re-scale
        propagator.parameters.Rs[1+3(i-2):3(i-1)] .= pos_3body
    end

    Uxx = propagator.parameters.f_jacobian(
        u[1:3]...,
        propagator.parameters.mus_scaled...,
        propagator.parameters.Rs...)
    return vcat(hcat(zeros(3,3), I(3)), hcat(reshape(Uxx, (3,3))', zeros(3,3)))
end


function jacobian(
    parameters::FullEphemParameters,
    et0::Float64,
    t::Real,
    u::Vector,
)
    # get third-body positions
    for i = 2:length(parameters.mus_scaled)
        # get position of third body
        pos_3body, _ = spkpos(
            parameters.naif_ids[i],
            et0 + t * parameters.tstar,
            parameters.naif_frame,
            parameters.abcorr,
            parameters.naif_ids[1]
        )
        pos_3body /= parameters.lstar                       # re-scale
        parameters.Rs[1+3(i-2):3(i-1)] .= pos_3body
    end

    Uxx = parameters.f_jacobian(
        u[1:3]...,
        parameters.mus_scaled...,
        parameters.Rs...)
    return vcat(hcat(zeros(3,3), I(3)), hcat(reshape(Uxx, (3,3))', zeros(3,3)))
end