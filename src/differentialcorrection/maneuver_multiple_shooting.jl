"""
Maneuver design multiple shooting
"""


"""
Forward-only multiple shooting problem
"""
mutable struct ManeuverForwardMultipleShootingProblem #<: AbstractDifferentialCorrectionProblem
    propagator#::PropagatorSTM
    epochs::Vector{Float64}
    nodes::Vector
    free_final_v::Bool
    maneuvers#::Vector{Float64}

    # constructor
    function ManeuverForwardMultipleShootingProblem(
        propagator,#::PropagatorSTM,
        epochs::Vector{Float64},
        nodes::Vector,
        free_final_v::Bool = true,
    )
        @assert length(epochs) == length(nodes) "Length of epochs and nodes should be equal!"

        new(
            propagator,
            epochs,
            nodes,
            free_final_v,
            Float64[]
        )
    end
end


"""
Solve multiple shooting problem.
The assumed variables are the nodes only.
"""
function shoot_fixedtime(problem::ManeuverForwardMultipleShootingProblem,
               ftol::Real=1e-8;
               maxiter::Int=1,
               verbose::Bool=true
    )

    # initialize storage
    N_nodes = length(problem.epochs)
    sols = Vector{ODESolution}(undef, N_nodes-1)
    DF = zeros(3(N_nodes-1), 3(N_nodes-1))
    residuals = zeros(3 * (N_nodes - 1))
    initial_node = deepcopy(problem.nodes[1])           # keep copy for later

    for idx in 1:maxiter
        for (idx, (et0, node)) in enumerate(zip(problem.epochs[1:end-1], problem.nodes[1:end-1]))
            # propagate each node forward in time
            tspan = (0.0, (problem.epochs[idx+1] - et0) / problem.propagator.parameters.tstar)
            sol = propagate(problem.propagator, et0, tspan, node)
            
            # position residuals
            residuals[1+3*(idx-1):3*idx] = problem.nodes[idx+1][1:3] - sol.u[end][1:3]
            sols[idx] = sol
            
            # Jacobian: ∂r(tf) / ∂v(t0)
            DF[1+3*(idx-1):3*idx, 1+3*(idx-1):3*idx] = -reshape(sol.u[end][7:42], (6,6))'[1:3,4:6]
        end
        
        if verbose
            @printf("Iteration %d: ||F|| = %e\n", idx, norm(residuals))
        end

        # check condition
        if norm(residuals) < ftol
            @printf("Achieved ||F|| = %e <= %e in %d iterations.\n", norm(residuals), ftol, idx)
            break
        end

        # compute Jacobian and perform update
        Δx = inv(DF) * residuals  #transpose(DF) * inv(DF * transpose(DF)) * residuals
        
        for idx in 1:N_nodes - 1
            problem.nodes[idx][4:6] -= Δx[1+3*(idx-1):3*idx]
        end
    end

    # store maneuvers
    problem.maneuvers = [[0.0, 0.0, 0.0] for _ in 1:N_nodes-1]  # initialize
    problem.maneuvers[1] .= sols[1].u[1][4:6] - initial_node[4:6]
    for (idx, (sol1, sol2)) in enumerate(zip(sols[1:end-1], sols[2:end]))
        problem.maneuvers[idx+1] .= sol2.u[1][4:6] - sol1.u[end][4:6]
    end
    return sols, residuals, DF
end