# Bivariate cooling crystallization with aggregation

This is the `(L,z)` bivariate migration of
`tutorials/oneWayCoupledCrysPBEFoam/shearIlievski`.  It retains the one-cell
K₂SO₄ cooling protocol, size-dependent growth, area-inhibited nucleation, and
aggregation.  The original crystallization-specific `crystalAggregation` /
`Ilievski` kernel is not part of the new generic bivariate kernel interface,
so this portable tutorial uses its existing `constant` aggregation kernel.

The cooling ramp drives supersaturation; competitive adsorption transfers
impurity from `Ci` into the second coordinate `z`.  The six transported
moments are `(0,0)`, `(1,0)`, `(2,0)`, `(3,0)`, `(0,1)`, and `(1,1)`.

Run `./Allrun` from an OpenFOAM environment containing `crysImpurityPBEFoam`.
