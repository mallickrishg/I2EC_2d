# Examples

The scripts in `examples/` run complete earthquake-cycle models on the example
subduction zone in `examples/inputs/`. Each runs from any folder, for example

```matlab
run examples/example1_numericalPowerLaw
```

and adds the I2EC folders to the path itself. All of them use the default
stress kernels: averaged over each mantle cell and certified free of growing
modes (see [Kernel options](../inputs.md#kernel-options) and
{ref}`mantle-kernels`), loaded from `examples/kernels/evl_orig.h5`.

| Example | Rheology | Solver | Run time |
| --- | --- | --- | --- |
| {doc}`example1` | power-law mantle ($n = 3$), rate-strengthening fault | numerical (`ode45`), 15-cycle spin-up | 1–2 min |
| {doc}`example2` | Newtonian mantle ($n = 1$), linear viscous fault | semi-analytical | 1–2 min |
| {doc}`example3` | as example 2, on three mantle meshes | semi-analytical | about 1 min (2 min the first time) |
| {doc}`example4` | as example 2 | both, compared | 1–2 min |
| {doc}`example5` | Maxwell, Burgers and Andrade mantles | semi-analytical | about 10 min |
| {doc}`example6` | Burgers mantle and fault inverted from surface displacement with `fmincon` | semi-analytical, adjoint | about 16 min |

The figures in these pages are made by `docs/figures/make_examples.m`, which
runs the scripts and saves what they draw.

```{toctree}
:maxdepth: 1
:hidden:

example1
example2
example3
example4
example5
example6
```

## The model

All examples share one geometry and loading: a megathrust dipping at 14°,
locked from the surface to 30 km depth and creeping below, a triangulated
oceanic mantle and mantle wedge (521 cells), and an earthquake every 500 years
that releases the slip deficit of the locked zone. The mantle is loaded by the
corner-flow strain rates of {mat:func}`computeLongTermStrainRates`, and the fault
by the plate rate.

To set up your own model step by step, see {doc}`../rheology`. Two helper
functions in `examples/` set up the example model:

- `exampleModel` loads the earth model, the megathrust, the model boundary,
  the kinematic hinge fault and a mantle mesh (`'ShearZone'`, default
  `examples/inputs/shearzone`). With `'StressKernels', true` or
  `'DisplacementKernels', true` it also returns the kernels: loaded from
  `examples/kernels/` when they match the meshes, computed otherwise.
- `newtonianCycle` sets the rheology of examples 2–4 (uniform viscosity
  $10^{19}$ Pa s, linear viscous fault with $\eta = 5\times10^{17}$ Pa s over
  the fault width, `rcv.faultLaw = 'linear'`) and solves the periodic cycle
  with {mat:func}`computePeriodicInitialCondition`. It returns the
  operator $A$, the forcing $c$, the periodic state, and a function
  `cyc.solve(tvec)` that evaluates slip rate, strain rates and their integrals
  at any times after the earthquake with {mat:func}`solveAnalyticIvp`.

## Long-term velocity

`plotLongtermVelocity` draws the long-term (plate-rate) velocity field that
loads the cycle, from {mat:func}`computeLongTermVelocity`:

```matlab
run examples/plotLongtermVelocity
```

```{figure} ../_static/examples/plotLongtermVelocity_1.png
:alt: Long-term velocity field of the example model
:width: 100%

Long-term velocity (arrows) over the whole model box: the oceanic plate
subducts at the plate rate and drives corner flow in the mantle wedge. The
mantle mesh is shaded, the megathrust is the thick line.
```
