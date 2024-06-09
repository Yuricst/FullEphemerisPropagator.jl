"""
Functions for conversions
"""


function LU2km(propagator::FullEphemPropagator, r)
    return r * propagator.parameters.lstar
end

function km2LU(propagator::FullEphemPropagator, r)
    return r / propagator.parameters.lstar
end

function VU2kms(propagator::FullEphemPropagator, v)
    return v * propagator.parameters.vstar
end

function kms2VU(propagator::FullEphemPropagator, v)
    return v / propagator.parameters.vstar
end

function TU2sec(propagator::FullEphemPropagator, t)
    return t * propagator.parameters.tstar
end

function sec2TU(propagator::FullEphemPropagator, t)
    return t / propagator.parameters.tstar
end

function dim2nondim(propagator, state::Vector{Float64})
    return state ./ [propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar]
end

function dim2nondim(propagator, state::SVector{6,Float64})
    return SA[
        state[1] / propagator.parameters.lstar,
        state[2] / propagator.parameters.lstar,
        state[3] / propagator.parameters.lstar,
        state[4] / propagator.parameters.vstar,
        state[5] / propagator.parameters.vstar,
        state[6] / propagator.parameters.vstar
    ]
end

function nondim2dim(propagator, state::Vector{Float64})
    return state .* [propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar]
end

function nondim2dim(propagator, state::SVector{6,Float64})
    return SA[
        state[1] * propagator.parameters.lstar,
        state[2] * propagator.parameters.lstar,
        state[3] * propagator.parameters.lstar,
        state[4] * propagator.parameters.vstar,
        state[5] * propagator.parameters.vstar,
        state[6] * propagator.parameters.vstar
    ]
end


"""Convert non-dimensional STM to dimensional STM"""
function nondim2dim_stm(propagator, stm)
    _stm = copy(stm)
    _stm[1:3, 4:6] *= propagator.parameters.lstar/propagator.parameters.vstar 
    _stm[4:6, 1:3] *= propagator.parameters.vstar/propagator.parameters.lstar 
    return _stm
end


"""Convert dimensional STM to non-dimensional STM"""
function dim2nondim_stm(propagator, stm)
    _stm = copy(stm)
    _stm[1:3, 4:6] /= propagator.parameters.lstar/propagator.parameters.vstar 
    _stm[4:6, 1:3] /= propagator.parameters.vstar/propagator.parameters.lstar 
    return _stm
end