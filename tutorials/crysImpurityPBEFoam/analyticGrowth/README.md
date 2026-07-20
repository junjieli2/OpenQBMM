# Analytic bivariate crystal growth

This one-cell case exercises the `(L,z)` crystal/impurity population balance
without spatial transport, nucleation, aggregation, breakup, or species
coupling. Competitive adsorption is non-zero and constant, so the two initial
quadrature nodes follow

```text
L(t) = L0 + (1 - theta) G t
z(t) = z0 + B theta/3 [L(t)^3 - L0^3]
```

where `G = CgGrowth sigma`,
`theta = Ki Ci/(1 + Ki Ci + Ks C)`, and
`B = eta rhoi surfaceFactor/mRef`.

Run `./Allrun`. It invokes `./Allclean` first, so every invocation rebuilds the
mesh and recomputes the solution instead of accepting stale logs or time
directories. `./Allclean` can also be used independently to remove generated
case output.

The standard-library-only `validate.py` script compares all six final moments
with this analytic solution and also checks adsorption and concentration
bounds.

The bivariate moment equations and the coupled `C`/`Ci` equations use an
internal first-order Euler time derivative. This is required by the split
source update and is independent of the `ddtSchemes/default` choice used by
the flow and temperature equations on a static mesh. For a dynamic mesh,
`ddt(U)` must also be `Euler` so its relative mesh flux satisfies the same
discrete geometric conservation law.
