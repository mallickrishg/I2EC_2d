# Theory and methods

I2EC models the earthquake cycle as small departures from steady plate motion.
The fault and the mantle are discretised into patches and triangular cells,
elastic kernels link them, and the result is a system of ordinary differential
equations that is either integrated in time or, when it is linear, solved
exactly. For linear models the same solution gives the gradient of a data
misfit with respect to viscosity. The full derivation is in
Mallick & Sathiakumar (2026), listed on the {doc}`../contributors` page.

| Page | What it covers |
| --- | --- |
| {doc}`equations` | Long-term plate motion and corner flow; how fault slip and mantle flow change the stress through elastic kernels |
| {doc}`rheology` | Power-law, Burgers and Andrade mantle; rate-strengthening, linear and Burgers fault; locking and earthquakes |
| {doc}`solvers` | Time integration with `ode45`, and the exact periodic solution of linear models by eigendecomposition |
| {doc}`kernels` | Half-space Green's functions, the boundary correction, and cell-averaged, certified mantle kernels |
| {doc}`gradients` | Adjoint gradients and tangents of a displacement misfit with respect to viscosity, for inversions |
| {doc}`limits` | Assumptions, the range of Maxwell times each solver handles, and practical guidance |

Read them in order: the loading and the equations, the rheology that closes
them, the solvers, then how the kernels are computed. The gradients build on
the semi-analytical solver.

```{toctree}
:maxdepth: 1
:hidden:

equations
rheology
solvers
kernels
gradients
limits
```
