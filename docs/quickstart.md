# Quick start

## Requirements

- MATLAB R2023b or later.
- The [MATLAB PDE Toolbox](https://www.mathworks.com/products/pde.html), only to
  build new meshes.
- The Statistics and Machine Learning Toolbox, only for {mat:func}`drchrnd`
  (`gamrnd`).

Clone the repository and run MATLAB from its root folder:

```bash
git clone https://github.com/mallickrishg/I2EC_2d.git
```

The driver scripts add the helper folder to the path and import the geometry
package themselves:

```matlab
addpath functions/
import('geometry.*')
```

Units are SI-like throughout: lengths in metres, time in seconds, stress in MPa.

## Examples

The `examples/` folder holds four examples of the same subduction zone: a
megathrust locked above 30 km depth with an earthquake every 500 years that
releases the slip deficit of the locked zone, a creeping fault below it, and a
viscous mantle loaded by corner flow. Run them from any folder:

```matlab
run examples/example2_analyticNewtonian
```

| Example | Solver | Rheology | Shows |
| --- | --- | --- | --- |
| `example1_numericalPowerLaw` | numerical (`ode45`), 15 spin-up cycles | power law, $n = 3$; rate-strengthening fault | slip rate through the cycle, strain-rate snapshots, surface velocity |
| `example2_analyticNewtonian` | semi-analytical, periodic cycle solved directly | Newtonian, $\eta = 10^{19}$ Pa s; linear viscous fault | strain-rate snapshots, slip rate, surface velocity |
| `example3_meshConvergence` | semi-analytical | as example 2 | strain-rate snapshots and profiles on three meshes (127, 495, 1955 cells) |
| `example4_analyticVsNumerical` | both | as example 2 | the two solvers agree to about $10^{-9}$ over the cycle |

`plotLongtermVelocity` plots the long-term (corner-flow) velocity field. The
shared helpers are `exampleModel` (paths, meshes, cached kernels) and
`newtonianCycle` (the model of examples 2-4). {doc}`examples/index` describes each
example in detail, with its figures.

Examples 1, 2 and 4 load the cached kernels in `kernels/` and take one to two
minutes. Example 3 builds and caches the kernels of its three meshes in
`kernels/meshconvergence/` on its first run (about 15 minutes, mostly the fine
mesh); later runs take about a minute.

### Semi-analytical solver (linear viscoelasticity)

For $n = 1$ and a linear viscous fault the cycle is a linear ODE,
$dy/dt = A\,y + c$, with $y$ the slip rate on the creeping patches and the
mantle strain rates. `examples/newtonianCycle.m` sets up example 2's model:

1. **Rheology**: a linear viscous fault (`rcv.Asigma` is the viscosity per unit
   length, `rcv.faultLaw = 'linear'`), mantle viscosity through `shz.alpha`,
   and the locked zone `rcv.pinnedPosition`.
2. **Long-term loading**: corner-flow strain rates from
   {mat:func}`getStrainratesLongterm`.
3. **Coseismic stress change** from slip on the locked zone.
4. **Operator**: $A = \eta^{-1} K$ with {mat:func}`assembleStressKernel`.
   With the default (certified) kernels $A$ has no growing modes and is used
   as it is.
5. **Periodic cycle**: the steady periodic initial condition from
   {mat:func}`computePeriodicInitialConditionBalanced`, evaluated at any time
   with {mat:func}`analyticSolveIVP`.

The core of the solve:

```matlab
A  = assembleStressKernel(evl, locked) ./ eta;      % no growing modes
y0 = computePeriodicInitialConditionBalanced(A, c, b, Trecur);
[Y, U] = analyticSolveIVP(A, c, y0, tvec);          % state and its time integral
[V, e22dot, e23dot] = extractComponentsSolutionVector(Y, locked, rcv.N, shz.N);
```

:::{important}
The default kernels are cell-averaged and certified, so the operator has no
growing modes and nothing is modified. {mat:func}`analyticSolveIVP` and
{mat:func}`computePeriodicInitialConditionBalanced` raise an error if the
operator has growing modes, which happens only with centroid kernels
(`'cellaverage', 1`). For those, negate them explicitly with
{mat:func}`negateGrowingModes` and build every forcing term that contains the
operator, such as the long-term loading `edot_lt = -A * y_pl`, with the negated
operator. See {ref}`caveats-keeping-the-system-damped`.
:::

### Numerical integration (general rheology)

Example 1 integrates the equations in time, which works for any rheology:

```matlab
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, stress_change, Ncycles, Trecur);
```

{mat:func}`runImposedEarthquakeCycles` starts from the long-term state and
repeats the imposed events `Ncycles` times, using {mat:func}`ode_viscoelastic`
and MATLAB's adaptive Runge-Kutta 4(5) solver `ode45` in
{mat:func}`viscoeq_sequence_simulation`; the last cycle is returned. Check that
it is periodic: example 1 prints how much the strain and slip rates change over
the last cycle. The fault law is set by `rcv.faultLaw`: `'ratestrengthening'`
(default, $d \log v/dt = \dot\tau/(a-b)\sigma$) or `'linear'` (linear viscous,
$\tau = \eta' v$, the law of the semi-analytical solver; example 4 uses it to
compare the two solvers). The default (certified) kernels are stable for any
power-law exponent, so no `kernelmodify` is needed.

(kernel-options)=
## Kernel options

{mat:func}`computeAllStressKernelsBem` builds certified kernels by default. The
earlier centroid kernels and their fixes for spurious growing modes remain
available:

| Call | Mantle stress | Growing modes | Use |
| --- | --- | --- | --- |
| `computeAllStressKernelsBem(rcv, shz, bnd)` | averaged over each cell, 3-point rule; certified | none, for any viscosity and power-law exponent | default, both solvers |
| `..., 'cellaverage', 7` | 7-point rule; certified | none | finer quadrature, about 2x the cost of 3-point |
| `..., 'cellaverage', 1, 'kernelcertify', false` | at the cell centroid | a few (43 on the example mesh) | earlier kernels; need one of the fixes below |
| `..., 'cellaverage', 1, 'kernelcertify', false, 'kernelmodify', 1` | at the centroid | zeroed (neutral) | earlier default for `ode45` runs |
| `..., 'cellaverage', 1, 'kernelcertify', false, 'kernelmodify', 2` | at the centroid | damped to $-\max\lvert\lambda\rvert$ | `ode45` runs with long spin-ups |

For the semi-analytical solvers with centroid kernels, negate the growing
modes of the operator explicitly: `[~, ~, A] = negateGrowingModes(A)`.
{mat:func}`modifyGrowingModes` applies `kernelmodify` 1 or 2 to an existing
kernel struct. See {ref}`caveats-keeping-the-system-damped` for why centroid
kernels have growing modes and how the treatments differ.

## Tests

From the repository root:

```matlab
addpath tests
runTests            % fast tests, under a minute
runTests('all')     % also rebuilds the example kernels and runs the examples
```

The suite checks the Green's functions, kernels, corner flow and solvers
against physical identities (free surface, dislocation jump, reciprocity,
Stokes equilibrium, the matrix exponential). Tests tagged `KnownIssue`
document bugs listed in `tasks.md` and are expected to fail until they are
fixed. See `tests/README.md` for how to read the results.

## Repository layout

| Path | Contents |
| --- | --- |
| `examples/` | Four examples (numerical $n = 3$, semi-analytical $n = 1$, mesh convergence, the two solvers compared) and the long-term velocity plot. |
| `functions/` | Solvers, kernel assembly, long-term loading and plotting helpers. |
| `+geometry/` | MATLAB package with the earth model and the fault and shear-zone classes. |
| `meshing/` | `createShearZoneMesh.m`: builds the triangular shear-zone mesh (PDE Toolbox). |
| `inputs/` | Example fault, boundary and shear-zone geometry files. |
| `kernels/` | Pre-computed stress (`evl_orig.mat`) and displacement (`devl.mat`) kernels. |
| `tests/` | Test suite (`matlab.unittest`); see `tests/README.md`. |
