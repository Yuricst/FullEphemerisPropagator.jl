module FullEphemerisPropagator

    using LinearAlgebra
    using SPICE
    using DifferentialEquations

    include("eoms.jl")
    include("propagator.jl")

    export Nbody_params, eom_Nbody_SPICE!
    export Propagator, propagate


end # module FullEphemerisPropagator
