# I2EC: Integral Imposed Earthquake Cycles (2D)

**I2EC** is a MATLAB code for viscoelastic earthquake cycles in
two-dimensional subduction zones. It models fault slip and power-law mantle
flow through repeated imposed earthquakes, and solves linear models
semi-analytically, orders of magnitude faster than time integration.

```{figure} _static/overview.png
:alt: Cross-section of the example model geometry
:width: 100%
:class: landing-figure

The example model: a megathrust locked to 30 km depth above a triangulated
oceanic mantle and mantle wedge.
```

:::{container} landing-cards

- **{doc}`Quick start <quickstart>`**

  Install, run the semi-analytical and numerical solvers, choose kernels.

- **{doc}`Set up the rheology <rheology>`**

  Put Maxwell, power-law, Burgers or Andrade rheology on the mesh and call a solver.

- **{doc}`Examples <examples/index>`**

  Five worked models with figures, from power-law flow and transient creep to mesh convergence.

- **{doc}`Theory and methods <theory>`**

  Governing equations, rheology, solvers and how the kernels are computed.

- **{doc}`Model inputs <inputs>`**

  Fault, boundary and mantle-mesh files, and the cached kernels.

- **{doc}`API reference <api/solvers>`**

  Every function, grouped by solver, loading, kernels and geometry.

- **{doc}`Cite I2EC <contributors>`**

  Authors, the paper to cite, and third-party components.

:::

```{toctree}
:caption: Guide
:maxdepth: 2
:hidden:

quickstart
rheology
examples/index
theory
inputs
```

```{toctree}
:caption: API reference
:maxdepth: 1
:hidden:

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
:hidden:

contributors
migration
GitHub repository <https://github.com/mallickrishg/I2EC_2d>
```
