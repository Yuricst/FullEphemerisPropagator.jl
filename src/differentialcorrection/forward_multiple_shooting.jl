"""
Forward multiple shooting differential correction
"""


"""
Convert list of epochs in seconds to list of intermediate times of flight, in TU
"""
function epochs2tofs(epochs, tstar)
    return [epochs[i+1] - epochs[i] for i in 1:length(epochs)-1] / tstar
end

"""
Convert list of times of flight in TU to list of epochs in seconds
"""
function tofs2epochs(et0, tofs, tstar)
    return vcat([et0], et0 .+ cumsum(tofs) * tstar)
end


"""
Forward-only multiple shooting problem
"""
mutable struct ForwardMultipleShootingProblem <: AbstractDifferentialCorrectionProblem
    propagator_state_only::Propagator
    propagator::PropagatorSTM
    epochs::Vector{Float64}
    nodes::Vector
    tofs::Vector

    # constructor
    function ForwardMultipleShootingProblem(
        propagator_state_only::Propagator,
        propagator::PropagatorSTM,
        epochs::Vector{Float64},
        nodes::Vector,
    )
        @assert length(epochs) == length(nodes) "Length of epochs and nodes should be equal!"
        new(
            propagator_state_only,
            propagator,
            epochs,
            nodes,
            epochs2tofs(epochs, propagator.parameters.tstar)
        )
    end
end


"""
Solve multiple shooting problem with fixed times.
The assumed variables are the nodes only.
"""
function shoot_fixedtime(problem::ForwardMultipleShootingProblem,
               ftol::Real=1e-8;
               maxiter::Int=1,
               verbose::Bool=true
    )

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
            DF[1+6*(idx-1):6*idx, 1+6*(idx-1):6*idx] = -reshape(sol.u[end][7:42], (6,6))'
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
        Δx = transpose(DF) * inv(DF * transpose(DF)) * residuals
        for idx in 1:N_nodes
            problem.nodes[idx] -= Δx[1+6*(idx-1):6*idx]
        end
    end
    return sols, residuals, DF
end



"""
Solve multiple shooting problem with free times.
The assumed variables are the nodes only.
"""
function shoot_freetime(problem::ForwardMultipleShootingProblem,
               ftol::Real=1e-8;
               maxiter::Int=1,
               verbose::Bool=true
    )

    # initialize storage
    N_nodes = length(problem.epochs)
    sols = Vector{ODESolution}(undef, N_nodes-1)
    DF = zeros(6(N_nodes-1), 6N_nodes + N_nodes - 1)
    residuals = zeros(6 * (N_nodes - 1))

    # set diagonals of DF to 1
    for i = 1:6(N_nodes - 1)
        DF[i,6 + i] = 1.0
    end
    
    for idx in 1:maxiter
        for (idx, (et0, node, tof)) in enumerate(zip(problem.epochs[1:end-1],
                                                     problem.nodes[1:end-1],
                                                     problem.tofs))
            # propagate each node forward in time
            tspan = (0.0, tof)
            sol = propagate(problem.propagator, et0, tspan, node)

            # store sensitivity to nodes
            residuals[1+6*(idx-1):6*idx] = problem.nodes[idx+1][1:6] - sol.u[end][1:6]
            sols[idx] = sol
            DF[1+6*(idx-1):6*idx, 1+6*(idx-1):6*idx] = -reshape(sol.u[end][7:42], (6,6))'
        end

        # function to compute sensitivity to tofs
        function func_residuals_vs_tofs(_tofs)
            # get epochs based on tofs
            _epochs = tofs2epochs(problem.epochs[1], _tofs, problem.propagator.parameters.tstar)
            # propagate each node forward in time
            _residuals = zeros(6 * (N_nodes - 1))
            for (idx, (_et0, node, _tof)) in enumerate(zip(_epochs[1:end-1],
                                                         problem.nodes[1:end-1],
                                                         _tofs))
                tspan = (0.0, _tof)
                sol = propagate(problem.propagator, _et0, tspan, node)
                _residuals[1+6*(idx-1):6*idx] = problem.nodes[idx+1][1:6] - sol.u[end][1:6]
            end
            return _residuals
        end
        DF[:,6N_nodes+1:end] = FiniteDifferences.jacobian(
            FiniteDifferences.central_fdm(5, 1),
            func_residuals_vs_tofs,
            problem.tofs,
        )[1]
        
        if verbose
            @printf("Iteration %d: ||F|| = %e\n", idx, norm(residuals))
        end

        # check condition
        if norm(residuals) < ftol
            @printf("Achieved ||F|| = %e <= %e in %d iterations.", norm(residuals), ftol, idx)
            break
        end

        # compute Jacobian and perform updates
        Δx = transpose(DF) * inv(DF * transpose(DF)) * residuals

        # update states
        for idx in 1:N_nodes
            problem.nodes[idx] -= Δx[1+6*(idx-1):6*idx]
        end

        # update tofs
        for idx in 1:N_nodes - 1
            problem.tofs[idx] -= Δx[6*N_nodes + idx]
        end
        
        # overwrite epochs
        problem.epochs = tofs2epochs(problem.epochs[1], problem.tofs, problem.propagator.parameters.tstar)
    end
    return sols, residuals, DF
end