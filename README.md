# FullEphemerisPropagator.jl

This repository contains an easily portable implementation of the full-ephemeris spacecraft dynamics (position & velocity). 
Using this set of code requires `SPICE.jl`, `DifferentialEquations.jl`,  `Symbolics`,  `SymbolicUtils`, and `LinearAlgebra.jl`.

- [DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/)
- [Step size control](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/)

Solvers to consider
- Taylor methods (https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#TaylorIntegration.jl)

## Quick example

For the N-body problem, we can first do some setup:

```julia
using SPICE
using DifferentialEquations

include("../src/FullEphemerisPropagator.jl")

# furnish spice kernels
spice_dir = ENV["SPICE"]   # modify as necessary

# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# define parameters
mus = [
    4.9028000661637961E+03,
    3.9860043543609598E+05,
    1.3271244004193938E+11,
]                                   # GMs
naif_ids = ["301", "399", "10"]     # NAIF IDs of bodies
naif_frame = "J2000"                # NAIF frame
abcorr = "NONE"                     # aberration  correction
lstar = 3000.0                      # canonical length scale
```

Now for integrating, there are two APIs available; the high-level API is as follows:

```julia
# instantiate propagator
prop = FullEphemerisPropagator.Propagator(
    Vern9(),
    lstar,
    mus,
    naif_ids;
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
u0 = [
    -2.5019204591096096,
    14.709398066624694,
    -18.59744250295792,
    5.62688812721852e-2,
    1.439926311669468e-2,
    3.808273517470642e-3
]

# time span (1 day, in canonical scale)
tspan = (0.0, 86400/prop.parameters.tstar)

# solve
sol = FullEphemerisPropagator.propagate(prop, et0, u0, tspan)
```

If it is desirable to use `DifferentialEquations.jl`'s calls to `ODEProblem()` and `solve()` directly, we can do:

```julia
# construct parameters
parameters = FullEphemerisPropagator.Nbody_params(
    et0,
    lstar,
    mus,
    naif_ids;
    naif_frame=naif_frame,
    abcorr=abcorr
)

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
u0 = [
    -2.5019204591096096,
    14.709398066624694,
    -18.59744250295792,
    5.62688812721852e-2,
    1.439926311669468e-2,
    3.808273517470642e-3
]

# time span (1 day, in canonical scale)
tspan = (0.0, 86400/parameters.tstar)

# solve
prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, parameters)
sol = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12)
```

Finally, plotting: 

```julia
using GLMakie
fig = Figure(resolution=(600,600), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=(1,1,1))
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])
fig
```


## To-do's

- [x] Jacobian via `Symbolics.jl`
- [ ] Spherical harmonics
- [ ] SRP
