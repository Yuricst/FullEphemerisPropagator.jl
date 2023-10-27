# FullEphemerisPropagator.jl

This repository contains an easily portable implementation of the full-ephemeris spacecraft dynamics (position & velocity). 
Using this set of code requires `SPICE.jl`, `DifferentialEquations.jl`, and `LinearAlgebra.jl`.

## Quick example

For the N-body problem, we can integrate as follows:

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

# initial epoch
et0 = str2et("2020-01-01T00:00:00")

# initial state (in canonical scale)
u0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

# time span (in canonical scale)
tspan = (0.0, 2*pi)
```

**High-level API**

```julia
# instantiate propagator
prop = FullEphemerisPropagator.Propagator(
    Tsit5(),
    lstar,
    mus,
    naif_ids;
    naif_frame = naif_frame,
    reltol = 1e-12,
    abstol = 1e-12,
)

# solve
sol = FullEphemerisPropagator.propagate(prop, et0, u0, tspan)
```

**Low-level API**

```julia
# construct parameters
params = FullEphemerisPropagator.Nbody_params(
    et0,
    lstar,
    mus,
    naif_ids;
    naif_frame=naif_frame,
    abcorr=abcorr
)

# solve
prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, params)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
```

Plotting: 

```julia
using GLMakie
fig = Figure(resolution=(600,600), fontsize=22)
ax1 = Axis3(fig[1, 1], aspect=(1,1,1))
lines!(ax1, sol[1,:], sol[2,:], sol[3,:])
fig
```


## To-do's

- [ ] Jacobian via `Symbolics.jl`
- [ ] Spherical harmonics
- [ ] SRP
