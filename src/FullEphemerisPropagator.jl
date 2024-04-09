module FullEphemerisPropagator

    using LinearAlgebra
    using StaticArrays
    using SPICE
    using OrdinaryDiffEq
    import Symbolics
    import SymbolicUtils
    import FiniteDifferences
    using Printf: @printf

    abstract type AbstractDifferentialCorrectionProblem end

    include("parameters.jl")
    include("eoms/perturbations.jl")
    include("eoms/eoms_Nbody.jl")
    include("eoms/eoms_NbodySRP.jl")
    include("symbolic_jacobians.jl")
    include("propagator.jl")
    include("canonical.jl")
    include("callbacks.jl")
    include("differentialcorrection/forward_multiple_shooting.jl")
    include("differentialcorrection/maneuver_multiple_shooting.jl")

    export Nbody_params, eom_Nbody_SPICE!, eom_Nbody_STM_SPICE!
    export Propagator, PropagatorSTM, propagate
    export LU2km, km2LU, VU2kms, kms2VU, TU2sec, sec2TU, dim2nondim, nondim2dim
    export ForwardMultipleShootingProblem, shoot_fixedtime, shoot_freetime

end
