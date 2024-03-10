# FullEphemerisPropagator.jl

`FullEphemerisPropagator.jl` is an astrodynamics library written in pure Julia for the propagation of the translational state of a spacecraft in high-fidelity dynamics models:

```math
\ddot{\boldsymbol{r}}
=
-\dfrac{\mu}{r^3} \boldsymbol{r}
+ \sum_i \boldsymbol{a}_{\mathrm{body}\,i}
+ \boldsymbol{a}_{\mathrm{srp}}
```

The library is essentially a wrapper around `DifferentialEquations.jl`, leveraging its numerical integration ecosystem. Ephemerides of celestial bodies are queried using JPL's SPICE kernels, which are handled via `SPICE.jl`. 

Current capabilities include:

- N-body equations of motion (restricted two-body + 3rd body perturbations)
- Solar radiation pressure
- Propagation of state-transition matrix (STM) together with the state

### Package dependencies

Using this package requires the following Julia packages:

- `LinearAlgebra`, `StaticArrays`, `SPICE`, `OrdinaryDiffEq`, `Symbolics`, `SymbolicsUtils`, `Printf`, `FiniteDifferences`

- Note: `FiniteDifferences` is used for performing differential correction, but the computation of STMs is done via Symbolic jacobians for speed & accuracy.
- In addition, relevant SPICE kernels must be downloaded and `furnsh`-ed. The minimum required kernels (for simple N-body propagation) is the leapseconds kernel (typically `naif0012.tls`, see [here](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/)) and the ephemeris kernel (typically `de440.bsp` for major solar system bodies, see [here](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/)). 
- Some generic kernels may be found at the NAIF website: [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/)

### Setup

To setup and use the library, run

```julia
pkg> add https://github.com/Yuricst/FullEphemerisPropagator.jl.git
```