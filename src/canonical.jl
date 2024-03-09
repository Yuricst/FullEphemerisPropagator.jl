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

function dim2nondim(propagator, state::Vector)
    return state ./ [propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar]
end

function nondim2dim(propagator, state::Vector)
    return state .* [propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.lstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar,
                     propagator.parameters.vstar]
end

