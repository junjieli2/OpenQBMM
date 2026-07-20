# Bivariate cooling crystallization with aggregation (3 L-nodes, 2 z-nodes)

This is the `(L,z)` bivariate crystallization case from `shearIlievski`
upgraded to a higher-order CQMOM closure with **3 L-nodes and 2 conditional
z-nodes** (`nNodes = [3, 2]`), using **15 transported moments**:

| Moment | Order (i, j) |
|--------|--------------|
| 1      | (0, 0)       |
| 2      | (1, 0)       |
| 3      | (2, 0)       |
| 4      | (3, 0)       |
| 5      | (4, 0)       |
| 6      | (5, 0)       |
| 7      | (0, 1)       |
| 8      | (1, 1)       |
| 9      | (2, 1)       |
| 10     | (0, 2)       |
| 11     | (1, 2)       |
| 12     | (2, 2)       |
| 13     | (0, 3)       |
| 14     | (1, 3)       |
| 15     | (2, 3)       |

The physics setup (K2SO4 cooling, size-dependent growth, area-inhibited
nucleation, constant aggregation kernel, competitive adsorption) is identical
to the `shearIlievski` case.  The higher-order closure provides improved
resolution of both the crystal length distribution and the conditional
impurity-content distribution.

Initial condition is monodisperse at L0 = 1e-4 m, z0 = 0.01.

Run `./Allrun` from an OpenFOAM environment containing `crysImpurityPBEFoam`.
