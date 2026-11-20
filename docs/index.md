# I2EC: Integral Imposed Earthquake Cycles (2D)

**I2EC** is a MATLAB code for viscoelastic earthquake cycles in two-dimensional
(plane-strain) subduction zones. An extended boundary element method turns the
governing equations into coupled ordinary differential equations for fault slip
rate and mantle strain rate, loaded by long-term plate motion and corner flow
and perturbed by periodic imposed earthquakes. The equations can be integrated
in time for general (power-law) rheology, or, for linear viscoelasticity,
solved semi-analytically by eigendecomposition, which is orders of magnitude
faster.

```{figure} _static/overview.png
:alt: Cross-section of the example model geometry
:width: 100%

The example model in `inputs/`: a megathrust dipping at 14° that is locked from
0 to 30 km depth and creeps below, above a triangulated oceanic mantle and below
a triangulated mantle wedge. The model boundary (`inputs/boundary2d.seg`) lies
outside this window. Made by `docs/figures/make_overview.m`.
```

```{toctree}
:caption: Guide
:maxdepth: 2

quickstart
examples/index
theory
inputs
```

```{toctree}
:caption: API reference
:maxdepth: 1

api/solvers
api/integration
api/longterm
api/kernels
api/geometry
api/utilities
```

```{toctree}
:caption: Project
:maxdepth: 1

contributors
GitHub repository <https://github.com/mallickrishg/I2EC_2d>
```
