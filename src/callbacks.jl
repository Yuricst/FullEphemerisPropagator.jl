"""
Helpful callback functions
"""

function get_maneuver_callbacks(
    et0::Float64,
    epochs::Vector{Float64},
    dv_vecs::Vector{Vector{Float64}},
    tstar::Float64,
    return_list::Bool = false,
)
    # define events for maneuvers
    _cbs = []
    for (et, dv_vec) in zip(epochs, dv_vecs)
        function _condition(u, t, integrator)
            return (et0 + t * tstar) - et
        end
        function _affect!(integrator)
            integrator.u[4:6] += dv_vec        # appending DV vector
            return 
        end
        push!(_cbs, ContinuousCallback(_condition, _affect!))
    end

    # return arguments
    if return_list
        return _cbs
    else
        return CallbackSet(_cbs...)
    end
end