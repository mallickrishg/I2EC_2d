# I2EC: Integral Imposed Earthquake Cycles (2D)

[![Docs](https://github.com/mallickrishg/I2EC_2d/actions/workflows/docs.yml/badge.svg)](https://mallickrishg.github.io/I2EC_2d/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs40623--026--02429--4-green)](https://doi.org/10.1186/s40623-026-02429-4)

I2EC is a MATLAB code for viscoelastic earthquake cycles in two-dimensional
(plane-strain) subduction zones. An extended boundary element method turns the
governing equations into coupled ordinary differential equations for fault slip
rate and mantle strain rate. These are loaded by long-term plate motion and
corner flow and perturbed by periodic imposed earthquakes.

**📖 Documentation: <https://mallickrishg.github.io/I2EC_2d/>**

The documentation is the reference for the theory, the input files, the
solvers and every function. This README only gets you started.

## Features

- **Two solvers.** Integrate in time with `ode45` for power-law mantle
  rheology and rate-strengthening faults, or, for linear viscoelasticity, solve
  the periodic cycle semi-analytically by eigendecomposition, which is orders
  of magnitude faster.
- **Consistent across timescales.** Cycle perturbations are added to long-term
  deformation (corner flow and steady fault slip), so the model stays
  consistent with plate motion over geological time.
- **Stable kernels by default.** Mantle stress kernels are averaged over each
  cell and certified free of spurious growing modes.
- **Predictions** of fault slip rate, mantle strain rate and stress, and
  surface displacement and velocity.

## Requirements

- MATLAB R2023b or later.
- [PDE Toolbox](https://www.mathworks.com/products/pde.html), only to build
  new meshes.
- Parallel Computing Toolbox is optional; kernel computation uses `parfor` if
  a pool is available.

## Quick start

```bash
git clone https://github.com/mallickrishg/I2EC_2d.git
```

In MATLAB, run any example from any folder (each adds the paths it needs):

```matlab
run I2EC_2d/examples/example2_analyticNewtonian
```

| Example | What it shows |
| --- | --- |
| [`example1_numericalPowerLaw`](https://mallickrishg.github.io/I2EC_2d/examples/example1.html) | Power-law ($n = 3$) mantle and rate-strengthening fault, integrated in time |
| [`example2_analyticNewtonian`](https://mallickrishg.github.io/I2EC_2d/examples/example2.html) | Newtonian mantle, periodic cycle solved semi-analytically |
| [`example3_meshConvergence`](https://mallickrishg.github.io/I2EC_2d/examples/example3.html) | The example 2 cycle on three mantle meshes |
| [`example4_analyticVsNumerical`](https://mallickrishg.github.io/I2EC_2d/examples/example4.html) | Both solvers on the same cycle |

Each takes 1–2 minutes (example 3 about 2 minutes the first time, while it
computes and caches its kernels). To build your own model, start from the
[quick start](https://mallickrishg.github.io/I2EC_2d/quickstart.html) and the
description of the [input files](https://mallickrishg.github.io/I2EC_2d/inputs.html).

## Repository layout

| Path | Contents |
| --- | --- |
| `functions/` | Solvers, kernel assembly, long-term loading, mesh generation, plotting helpers |
| `+geometry/` | Earth model, fault and shear-zone classes, Green's functions |
| `examples/` | The four examples and their shared model setup |
| `examples/inputs/` | Geometry of the example subduction zone |
| `examples/kernels/` | Pre-computed kernels of the example model (HDF5) |
| `tests/` | `matlab.unittest` test suite |
| `docs/` | Source of the documentation website |
| `*.nb` | Mathematica notebooks for the corner-flow solutions |

## Where to find what

| I want to… | Go to |
| --- | --- |
| Understand the equations and the method | [Theory](https://mallickrishg.github.io/I2EC_2d/theory.html) |
| Set up my own fault, mesh or kernels | [Inputs](https://mallickrishg.github.io/I2EC_2d/inputs.html) |
| Choose kernel options | [Quick start: kernel options](https://mallickrishg.github.io/I2EC_2d/quickstart.html#kernel-options) |
| Look up a function | [API reference](https://mallickrishg.github.io/I2EC_2d/api/solvers.html) |
| See worked examples with figures | [Examples](https://mallickrishg.github.io/I2EC_2d/examples/index.html) |

## Tests

```matlab
addpath tests
runTests          % fast tests, about a minute
runTests('all')   % also rebuilds the example kernels and runs the examples
```

See [`tests/README.md`](tests/README.md) for how to read the results. Open
maintenance items are tracked in [`tasks.md`](tasks.md).

## Contributing

Bug reports and questions are welcome as
[GitHub issues](https://github.com/mallickrishg/I2EC_2d/issues). For code
changes, open a pull request against `main` and run `runTests` first.

## Citation

If you use I2EC in your research, please cite:

> Mallick, R., & Sathiakumar, S. (2026). A fast numerical routine to model
> viscoelastic earthquake cycles in subduction zones. *Earth, Planets and
> Space*. <https://doi.org/10.1186/s40623-026-02429-4>

## Authors and license

Rishav Mallick and Sharadha Sathiakumar (Earth Observatory of Singapore); see
[contributors](https://mallickrishg.github.io/I2EC_2d/contributors.html) for
who did what and the third-party components I2EC builds on.

Released under the [MIT License](LICENSE).
