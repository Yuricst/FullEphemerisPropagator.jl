module FullEphemerisPropagator

    using LinearAlgebra
    using SPICE
    using DifferentialEquations
    import Symbolics
    import SymbolicUtils
    using Printf: @printf

    include("eoms.jl")
    include("propagator.jl")
    include("differentialcorrection.jl")

    export Nbody_params, eom_Nbody_SPICE!, eom_NbodySTM_SPICE!
    export Propagator, PropagatorSTM, propagate
    export ForwardMultipleShootingProblem, shoot

end # module FullEphemerisPropagator
