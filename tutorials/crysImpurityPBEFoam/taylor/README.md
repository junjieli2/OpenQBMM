# Bivariate Taylor–Couette crystallization

This is the bivariate `(L,z)` counterpart of
`tutorials/oneWayCoupledCrysPBEFoam/taylor`.  It retains the rotating-cylinder
Taylor–Couette geometry and couples flow, heat transport, solute transport,
impurity transport, nucleation, growth, and competitive impurity adsorption.

The mesh and run duration are deliberately reduced so that the tutorial is a
quick smoke test.  The six moments describe crystal length `L` and incorporated
impurity content `z`:

```text
(0,0) (1,0) (2,0) (3,0) (0,1) (1,1)
```

Run `./Allrun` from an OpenFOAM environment containing `crysImpurityPBEFoam`.
