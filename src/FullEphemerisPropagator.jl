module FullEphemerisPropagator

    using LinearAlgebra
    using SPICE
    using DifferentialEquations

    include("eoms.jl")


    export Nbody_params, eom_Nbody_SPICE!


end # module FullEphemerisPropagator
