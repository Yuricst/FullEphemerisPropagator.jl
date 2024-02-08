module FullEphemerisPropagator

    using LinearAlgebra
    using SPICE
    using DifferentialEquations
    import Symbolics
    import SymbolicUtils

    #include("symbolics/f_jacobian_NbodyN2.jl")

    include("eoms.jl")
    include("propagator.jl")

    export Nbody_params, eom_Nbody_SPICE!, eom_NbodySTM_SPICE!
    export Propagator, PropagatorSTM, propagate


end # module FullEphemerisPropagator
