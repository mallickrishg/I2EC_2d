# Quick start

## Install

You need MATLAB R2023b or later. The
[PDE Toolbox](https://www.mathworks.com/products/pde.html) is needed only to
build new meshes.

```bash
git clone https://github.com/mallickrishg/I2EC_2d.git
```

In your own scripts, put the code on the path from the repository root:

```matlab
addpath functions/
import('geometry.*')
```

Units: lengths in metres, time in seconds, stress in MPa.

## Run your first model

```matlab
run examples/example2_analyticNewtonian
```

In about a minute this solves a full earthquake cycle of an example subduction
zone and plots the mantle strain rate, the slip rate on the fault and the
surface velocity at three times after the earthquake. {doc}`examples/example2`
walks through the script step by step.

## Choose an example

| Example | What it shows |
| --- | --- |
| {doc}`examples/example1` | Power-law mantle and friction, integrated in time |
| {doc}`examples/example2` | Newtonian mantle, cycle solved exactly (fastest) |
| {doc}`examples/example3` | How the answer depends on the mesh |
| {doc}`examples/example4` | The two solvers agree on the same model |
| {doc}`examples/example5` | Burgers and Andrade transient creep compared with Maxwell |
| {doc}`examples/example6` | Invert surface displacement for mantle viscosity (adjoint gradients) |

All six use the same geometry and cached kernels; see {doc}`examples/index`.

## Build your own model

1. **Geometry:** write the fault, boundary and mantle-mesh files
   ({doc}`inputs`).
2. **Kernels:** compute them once with
   {mat:func}`computeAllStressKernelsBem` and save them
   ({doc}`inputs`, [kernel options](inputs.md#kernel-options)).
3. **Rheology and solve:** set the fault and mantle properties and call a
   solver, step by step in {doc}`rheology`. Check that your mantle
   Maxwell time is in the range the solver handles ({ref}`valid-range`).

{doc}`theory/index` explains the equations behind each step.

## Check your installation

```matlab
addpath tests
runTests            % about a minute
```

See `tests/README.md` for how to read the results.
