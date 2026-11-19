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

## Semi-analytical solver (linear viscoelasticity)

`MAIN_analyticsolve_singleeqcycle.m` is the fastest way to get results. It loads
pre-computed kernels from `kernels/`, so no boundary-element assembly is needed:

```matlab
run MAIN_analyticsolve_singleeqcycle
```

The script steps through:

1. **Geometry**: the earth model {mat:class}`geometry.LDhs <LDhs>`, the
   megathrust and model boundary ({mat:class}`geometry.receiver <receiver>`) and
   the triangulated mantle
   ({mat:class}`geometry.shearZoneReceiver <shearZoneReceiver>`).
2. **Kernels**: loads `evl_orig` (stress) and `devl` (surface displacement). To
   recompute them, use {mat:func}`computeAllStressKernelsBem` and
   {mat:func}`computeAllDisplacementKernelsBem`.
3. **Rheology**: an effective fault viscosity `rcv.Asigma`, mantle viscosities
   through `shz.alpha`, and the locked zone `rcv.pinnedPosition`.
4. **Long-term loading**: corner-flow strain rates from
   {mat:func}`getStrainratesLongterm`.
5. **Coseismic stress change** from slip on the locked zone.
6. **Operator**: assembles $A = \eta^{-1} K$ with {mat:func}`assembleStressKernel`
   and negates spurious growing modes with {mat:func}`negateGrowingModes`.
7. **Periodic cycle**: the steady periodic initial condition from
   {mat:func}`computePeriodicInitialConditionBalanced`.
8. **Postseismic evolution** at chosen times with {mat:func}`analyticSolveIVP`,
   split into fault slip rate and strain rates with
   {mat:func}`extractComponentsSolutionVector`.
9. **Plots** of strain-rate snapshots, surface velocities and station time series.

The core of the solve:

```matlab
stresskernel = assembleStressKernel(evl_orig, locked);
rheoparam    = diag(1./viscosityvector) * stresskernel;
[~, ~, rheoparam] = negateGrowingModes(rheoparam);

c  = edot_locking + edot_lt;          % forcing
y0 = computePeriodicInitialConditionBalanced(rheoparam, c, deltastrainrate, Trecur);

tvec = [40, 365, 10*365, 50*365] .* 86400;
[sol, sol_integrated] = analyticSolveIVP(rheoparam, c, y0, tvec);
[V, e22dot, e23dot] = extractComponentsSolutionVector(sol, locked, rcv.N, shz.N);
```

:::{important}
The kernels and the operator are modified to keep the system damped; see
{ref}`caveats-keeping-the-system-damped` before interpreting cell-scale
features. Build every forcing term that contains the operator, such as the long-term
loading `edot_lt = -A * y_pl`, with the *negated* operator, so that the
plate-rate state stays the steady state of the evolved system. See
{mat:func}`negateGrowingModes`.
:::

## Numerical integration (general rheology)

`MAIN_i2ec_singleeqcycle.m` computes the kernels from scratch and integrates the
equations in time:

```matlab
run MAIN_i2ec_singleeqcycle
```

After setting up geometry, kernels, rheology and coseismic stress changes as
above, it calls:

```matlab
Ncycles = 5;   % cycles to spin up to a periodic state
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, stress_change, Ncycles, Trecur);
```

{mat:func}`runImposedEarthquakeCycles` repeats the imposed events until the
response is periodic, using {mat:func}`ode_viscoelastic` and MATLAB's adaptive
Runge-Kutta 4(5) solver `ode45` in {mat:func}`viscoeq_sequence_simulation`. Use this path
for power-law rheology, where the semi-analytical solver does not apply.

## Tests

From the repository root:

```matlab
addpath tests
runTests            % fast tests, under a minute
runTests('all')     % also rebuilds the example kernels and runs the driver scripts
```

The suite checks the Green's functions, kernels, corner flow and solvers
against physical identities (free surface, dislocation jump, reciprocity,
Stokes equilibrium, the matrix exponential). Tests tagged `KnownIssue`
document bugs listed in `tasks.md` and are expected to fail until they are
fixed. See `tests/README.md` for how to read the results.

## Repository layout

| Path | Contents |
| --- | --- |
| `MAIN_i2ec_singleeqcycle.m` | Driver: imposed earthquake cycles by numerical integration (`ode45`). |
| `MAIN_analyticsolve_singleeqcycle.m` | Driver: linear viscoelastic cycles with the semi-analytical solver. |
| `PLOT_longtermVelocity.m` | Plots the long-term (corner-flow) velocity field. |
| `functions/` | Solvers, kernel assembly, long-term loading and plotting helpers. |
| `+geometry/` | MATLAB package with the earth model and the fault and shear-zone classes. |
| `meshing/` | Scripts to build triangular shear-zone meshes (based on *mesh2d*). |
| `inputs/` | Example fault, boundary and shear-zone geometry files. |
| `kernels/` | Pre-computed stress (`evl_orig.mat`) and displacement (`devl.mat`) kernels. |
| `tests/` | Test suite (`matlab.unittest`); see `tests/README.md`. |
