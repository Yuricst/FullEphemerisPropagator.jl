# FullEphemerisPropagator.jl

This repository contains an easily portable implementation of the full-ephemeris spacecraft dynamics (position & velocity). 
Using this set of code requires `SPICE.jl`, `DifferentialEquations.jl`, and `LinearAlgebra.jl`.

## Quick example

For the N-body problem, we can integrate as follows:

```julia
using SPICE
using Plots
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
]
naif_ids = ["301", "399", "10"]
naif_frame = "J2000"
abcorr = "NONE"
lstar = 3000.0

et0 = str2et("2020-01-01T00:00:00")
params = FullEphemerisPropagator.Nbody_params(
    et0,
    lstar,
    mus,
    naif_ids;
    naif_frame=naif_frame,
    abcorr=abcorr
)

# initial state (in canonical scale)
u0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

# time span (in canonical scale)
tspan = (0.0, 2*pi)

# solve
prob = ODEProblem(FullEphemerisPropagator.eom_Nbody_SPICE!, u0, tspan, params)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)

# plot
plot(sol, idxs=(1,2,3), label="Trajectory")
```


## To-do's

- [ ] Jacobian via `Symbolics.jl`
- [ ] Spherical harmonics
- [ ] SRP
