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

