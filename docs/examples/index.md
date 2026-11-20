# Examples

The scripts in `examples/` run complete earthquake-cycle models on the example
subduction zone in `inputs/`. Each runs from any folder, for example

```matlab
run examples/example2_analyticNewtonian
```

and adds the I2EC folders to the path itself. All of them use the default
stress kernels: averaged over each mantle cell and certified free of growing
modes (see [Kernel options](quickstart.md#kernel-options) and
{ref}`caveats-keeping-the-system-damped`), loaded from `kernels/evl_orig.mat`.

| Script | Rheology | Solver | Run time |
| --- | --- | --- | --- |
| `example1_numericalPowerLaw` | power-law mantle ($n = 3$), rate-strengthening fault | numerical (`ode45`), 15-cycle spin-up | 1–2 min |
| `example2_analyticNewtonian` | Newtonian mantle ($n = 1$), linear viscous fault | semi-analytical | seconds |
| `example3_meshConvergence` | as example 2, on three mantle meshes | semi-analytical | about 1 min (about 15 min the first time) |
| `example4_analyticVsNumerical` | as example 2 | both, compared | 1–2 min |
| `plotLongtermVelocity` | corner flow | — | seconds |

The figures on this page are made by `docs/figures/make_examples.m`, which runs
the scripts and saves what they draw.

## The model

All examples share the geometry and loading of {doc}`quickstart`: a megathrust
dipping at 14°, locked from the surface to 30 km depth and creeping below, a
triangulated oceanic mantle and mantle wedge (521 cells), and an earthquake
every 500 years that releases the slip deficit of the locked zone. The
mantle is loaded by the corner-flow strain rates of
{mat:func}`getStrainratesLongterm`, and the fault by the plate rate.

Two helper functions set this up:

- `exampleModel` loads the earth model, the megathrust, the model boundary,
  the kinematic hinge fault and a mantle mesh (`'ShearZone'`, default
  `inputs/shearzone`). With `'StressKernels', true` or
  `'DisplacementKernels', true` it also returns the kernels: loaded from
  `kernels/` when they match the meshes, computed otherwise.
- `newtonianCycle` sets the rheology of examples 2–4 (uniform viscosity
  $10^{19}$ Pa s, linear viscous fault with $\eta = 5\times10^{17}$ Pa s over
  the fault width, `rcv.faultLaw = 'linear'`) and solves the periodic cycle
  with {mat:func}`computePeriodicInitialConditionBalanced`. It returns the
  operator $A$, the forcing $c$, the periodic state, and a function
  `cyc.solve(tvec)` that evaluates slip rate, strain rates and their integrals
  at any times after the earthquake with {mat:func}`analyticSolveIVP`.

`plotLongtermVelocity` draws the long-term (plate-rate) velocity field that
loads the cycle, from {mat:func}`computeLongTermVelocity`:

```{figure} _static/examples/plotLongtermVelocity_1.png
:alt: Long-term velocity field of the example model
:width: 100%

Long-term velocity field (arrows) over the whole model box, with the mantle
mesh (red), the megathrust and the model boundary.
```

(example-2)=
## Example 2: Newtonian mantle, semi-analytical solver

`example2_analyticNewtonian` is the simplest complete model. The mantle is
Newtonian and the fault linear viscous, so the cycle is a linear ODE
$dy/dt = Ay + c$. Its periodic state is found directly, with no spin-up, and
the solution is evaluated exactly at $t/T = 0.01$, $0.1$ and $0.99$ after the
earthquake. The script also prints the slowest and fastest relaxation times
of $A$.

```{figure} _static/examples/example2_analyticNewtonian_1.png
:alt: Mantle strain rate long-term and at three times after the earthquake
:width: 100%

Mantle strain rate: long-term (top left) and at three times after the
earthquake. Shortly after it the whole mantle deforms faster than in the long
term; the excess decays through the cycle.
```

```{figure} _static/examples/example2_analyticNewtonian_2.png
:alt: Slip rate on the megathrust at three times
:width: 80%

Slip rate on the creeping part of the megathrust (the locked zone, $v = 0$,
is not shown on the log scale).
```

```{figure} _static/examples/example2_analyticNewtonian_3.png
:alt: Horizontal and vertical surface velocity at three times
:width: 80%

Horizontal (top) and vertical (bottom) surface velocity, normalised by the
plate rate.
```

(example-1)=
## Example 1: power-law mantle, numerical solver

`example1_numericalPowerLaw` keeps the geometry, locking and earthquake of
example 2 but makes the mantle a single power-law material,
$\dot\varepsilon = \alpha\,\tau^{n-1} s$ with $n = 3$ and uniform $\alpha$,
and the creeping fault rate-strengthening,
$d\log v/dt = \dot\tau/(a-b)\sigma$ with $(a-b)\sigma = 0.2$ MPa. Because
$\alpha$ is uniform, the $n = 3$ corner flow is the long-term solution of this
mantle; $\alpha$ is chosen so that the effective viscosity is $10^{19}$ Pa s at
the median long-term strain rate. The coseismic stress step on the creeping
fault is capped at 5 MPa next to the locked zone.

The equations are nonlinear, so {mat:func}`runImposedEarthquakeCycles`
integrates them in time with `ode45`, from the long-term state through
15 cycles. The script prints the change over the last cycle as a check of
periodicity: at most $3.4\times10^{-5}$ per cell (median $2\times10^{-7}$).
The long-term effective viscosity then ranges from $3\times10^{17}$ to
$7\times10^{19}$ Pa s.

```{figure} _static/examples/example1_numericalPowerLaw_1.png
:alt: Slip rate on the megathrust over the last cycle
:width: 80%

Slip rate on the megathrust over the last cycle (log time axis). Below the
locked zone ($x > 120$ km) afterslip starts at up to 100 times the plate rate
and decays over the cycle.
```

```{figure} _static/examples/example1_numericalPowerLaw_2.png
:alt: Mantle strain rate long-term and at three times in the last cycle
:width: 100%

Mantle strain rate in the last cycle. With $n = 3$ the stressed regions
weaken, so the postseismic strain rates decay faster than in example 2: by
$t/T = 0.1$ they are close to their late-interseismic level.
```

```{figure} _static/examples/example1_numericalPowerLaw_3.png
:alt: Surface velocity at three times in the last cycle
:width: 80%

Surface velocity in the last cycle.
```

(example-3)=
## Example 3: mesh convergence

`example3_meshConvergence` solves the cycle of example 2 on three mantle
meshes of the same domain, made with `meshing/createShearZoneMesh` with
largest edges of 40, 20 and 10 km (127, 495 and 1955 cells). The fault and
boundary meshes are the same, and the semi-analytical solver is exact in
time, so any difference between the three solutions comes from the mantle
mesh. The script prints the relative difference from the fine mesh at
$t/T = 0.01$, $0.1$ and $0.99$ and plots strain-rate snapshots and profiles.

The kernels of the three meshes are cached in `kernels/meshconvergence/`
(git-ignored). The first run computes them, about 15 minutes, mostly for the
fine mesh; later runs take about a minute.

(example-4)=
## Example 4: semi-analytical against numerical

`example4_analyticVsNumerical` solves the cycle of example 2 with both
solvers. The numerical run uses the linear viscous fault law
(`rcv.faultLaw = 'linear'`, set by `newtonianCycle`), so both solve the same
linear ODE. It starts from the analytic state just before the earthquake and
integrates one cycle with {mat:func}`viscoeq_sequence_simulation` and
{mat:func}`ode_viscoelastic`. It must follow the analytic solution at every
time and end where it started. On the example they agree to
$3.8\times10^{-9}$ in slip rate and $3.6\times10^{-11}$ in strain rate
(relative), and the numerical end state matches its start to
$3\times10^{-11}$.

```{figure} _static/examples/example4_analyticVsNumerical_1.png
:alt: Time series of fault slip rate and mantle strain rate from both solvers
:width: 100%

Slip rate at three creeping patches (32, 40 and 48 km depth) and strain rate
at three mantle cells through one cycle: semi-analytical (coloured) and
numerical (black dashed).
```

Starting the numerical run from the long-term state instead would need a
spin-up. The slowest mode of this model decays over thousands of years
(printed by the script), so it takes many 500-year cycles to become
periodic. `TimeIntegrationTest/spunUp*` checks the spin-up on a smaller model.
