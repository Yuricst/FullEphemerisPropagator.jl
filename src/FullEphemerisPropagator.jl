module FullEphemerisPropagator

    using LinearAlgebra
    using SPICE
    using OrdinaryDiffEq
    import Symbolics
    import SymbolicUtils
    using Printf: @printf

    include("eoms.jl")
    include("symbolic_jacobians.jl")
    include("propagator.jl")
    include("canonical.jl")
    include("differentialcorrection.jl")

    export Nbody_params, eom_Nbody_SPICE!, eom_NbodySTM_SPICE!
    export Propagator, PropagatorSTM, propagate
    export LU2km, km2LU, VU2kms, kms2VU, TU2sec, sec2TU, dim2nondim, nondim2dim
    export ForwardMultipleShootingProblem, shoot

end
