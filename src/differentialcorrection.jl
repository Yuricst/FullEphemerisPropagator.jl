"""
Differential correction functions
"""

abstract type AbstractDifferentialCorrectionProblem end


"""
Forward-only multiple shooting problem
"""
mutable struct ForwardMultipleShootingProblem <: AbstractDifferentialCorrectionProblem
    propagator::PropagatorSTM
    epochs::Vector{Float64}
    nodes::Vector

    # constructor
    function ForwardMultipleShootingProblem(
        propagator::PropagatorSTM,
        epochs::Vector{Float64},
        nodes::Vector,
    )
        @assert length(epochs) == length(nodes) "Length of epochs and nodes should be equal!"
        new(
            propagator,
            epochs,
            nodes,
        )
    end
end


"""
Solve multiple shooting problem.
The assumed variables are the nodes only.
"""
function shoot(problem::ForwardMultipleShootingProblem, ftol::Real=1e-8;
               maxiter::Int=1, verbose::Bool=true)

    # initialize storage
    N_nodes = length(problem.epochs)
    sols = Vector{ODESolution}(undef, N_nodes-1)
    DF = zeros(6(N_nodes-1), 6N_nodes)
    residuals = zeros(6 * (N_nodes - 1))

    # set diagonals of DF to 1
    for i = 1:6(N_nodes - 1)
        DF[i,6 + i] = 1.0
    end
    
    for idx in 1:maxiter
        for (idx, (et0, node)) in enumerate(zip(problem.epochs[1:end-1], problem.nodes[1:end-1]))
            # propagate each node forward in time
            tspan = (0.0, (problem.epochs[idx+1] - et0) / problem.propagator.parameters.tstar)
            sol = propagate(problem.propagator, et0, tspan, node)
            # store
            residuals[1+6*(idx-1):6*idx] = problem.nodes[idx+1][1:6] - sol.u[end][1:6]
            sols[idx] = sol
            DF[1+6*(idx-1):6*idx, 1+6*(idx-1):6*idx] = -reshape(sol.u[end][7:42], (6,6))
        end
        
        if verbose
            @printf("Iteration %d: ||F|| = %e\n", idx, norm(residuals))
        end

        # check condition
        if norm(residuals) < ftol
            @printf("Achieved ||F|| = %e <= %e in %d iterations.", norm(residuals), ftol, idx)
            break
        end

        # compute Jacobian and perform update
        Δx = transpose(DF) * inv(DF * transpose(DF)) * residuals
        for idx in 1:N_nodes
            problem.nodes[idx] -= Δx[1+6*(idx-1):6*idx]
        end
    end
    return sols, residuals, DF
end