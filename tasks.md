# Maintenance tasks

Code-maintenance backlog for I2EC_2d. Findings come from a static reference scan of every tracked `.m` file, the notebook (`I2EC_linearVEcycles.ipynb`, since deleted) and the Sphinx API pages (`docs/api/*.rst`), done on 2026-11-19.

"Unused" below means **nothing in the repo calls it** (comments excluded). Before deleting anything, check that no external experiment scripts or collaborators rely on it.

When you remove or rename a function, also update the `mat:autofunction::` entry in `docs/api/*.rst`. Otherwise the docs build will warn or fail.

**Run the test suite before and after each change** (see [`tests/README.md`](tests/README.md)):

```matlab
addpath tests; runTests          % fast tier
addpath tests; runTests('all')   % + kernel rebuild and MAIN scripts
```

Tests tagged `KnownIssue` check the *correct* behaviour of the bugs listed in §0. They fail today. When one starts passing, remove its tag.

## Still open (2026-11-21, checked against `main` after PR #31)

Everything not listed here is done; details are in the sections named.

**Code**
- [ ] §6: a regeneration script, or a mesh-size/hash check, for `examples/kernels/evl_orig.h5` and `devl.h5`.
- [ ] §7: one `setup.m` that adds the I2EC folders to the path, called by the examples and the tests.
- [ ] §8: run `runTests` in CI, and run `docs/_ext/check_help.py` in the docs workflow.

**Optional**
- [ ] §0a: rename `evl_orig` (it is no longer an "original, unmodified" kernel).
- [ ] §3: shorten `filenameGenerator.m` to `sprintf('%06d', num)`.
- [ ] §11: the stress generator reads the displacement file by line number; regenerate after any edit to `computeDisplacementPlaneStrainTriangleShearZone`.

**Future**
- [ ] §0: a particular solution for a prescribed boundary motion (`bcx`/`bcz ≠ 0`).

**Features (in priority order, 2026-11-21)**

Each feature gets its own branch from `main`, its own worktree outside Dropbox, and one agent. Its plan and checklist live in `docs/dev/<feature>.md` on that branch; this list keeps only one line per feature. Features never merge into each other: each starts from (or merges) `main` once the features before it have landed. Do not check out, commit to, merge or rebase another feature's branch.

1. - [ ] **Burgers/Andrade transient rheology** as an option in both the semi-analytical solver (extra transient-strain state per cell; Andrade as a Prony series of Kelvin elements) and the numerical solver. *In progress:* session `burgers-andrade-rheology`, notes in `docs/dev/transient-rheology.md`.
2. - [ ] **Gravity as a restoring force (buoyancy)**, for Maxwell and Burgers/Andrade. Buoyancy depends on displacement, so the state gains accumulated displacement and the linear periodic solve still applies. It must act on departures from the long-term topography. *Starts after 1 is merged.*
3. - [ ] **Adjoint operator/ODE for inverse problems** (e.g. estimating η(x) from surface velocities), including the adjoint of the periodicity condition. It covers every rheology built by then, plus gravity. *Starts after 2 is merged.*
4. - [ ] **Nonlinear (power-law) spin-up without brute force** (§10): Newton shooting on the post-earthquake state. Lowest priority, *paused*. Lives only on branch `nonlin-spinup-solve` (plan committed, `main` merged in), which is reserved for one agent.

**Done since 2026-11-20** (details in the sections): deviatoric `LL`/`LK` construction verified (§1); example 3 figures and caches; HDF5 kernel files; analytic triangle stress (§11, #14/#15); coupled fault-mantle check, closed with no code change (§0a); both safeguards raise errors instead of fixing quietly, and the solvers eigendecompose without balancing (§0a); legacy bugs and unused functions deleted (§0, §3); one corner-flow implementation (§1); naming rules, renames and third-party credits (§2, §4, #30/#31); example inputs, kernels and the mesher moved under `examples/` and `functions/`; help-block format (§9).

---

## 0. Bugs confirmed by the test suite (2026-11-19)

Each item names the failing test (`Class/test`).

### Affects code that is in use

- [x] ~~**`MAIN_i2ec_singleeqcycle.m` and `computeInterseismicSteadystate.m` cannot run.** *Fixed: `n` defaults to 1, and both callers now pass `n = 1`.* Both call `getStrainratesLongterm(shz, dip, oc, co)` with 4 arguments, but the function now requires the exponent `n`, so they fail with "Not enough input arguments". → Give `n` a default of 1 in `getStrainratesLongterm`, or pass it at both call sites. Tests: `LongTermFlowTest/strainRatesWithoutExponentDefaultToNewtonian`, `SolversTest/interseismicSteadyStateRuns`, `ExampleScriptsTest/numericalSingleCycleScriptRuns`.~~
- [x] ~~**`runImposedEarthquakeCycles(..., nreps = 1, ...)` fails.** The first call discards `t` (`[\~,Y] = ...`) and the loop never runs, so output `t` is never assigned. → Keep `t` from the first call. *Fixed: the first call now keeps `t`, so `nreps = 1` returns one cycle from the long-term state. `nreps > 1` is unchanged. Test: `TimeIntegrationTest/singleCycleRunReturnsTime` (no longer `KnownIssue`; also checks `t(end) = Trecur` and the output sizes). Since then `viscoeq_sequence_simulation` takes `nreps` and runs every segment of every cycle in one loop, so `runImposedEarthquakeCycles` no longer has a cycle loop of its own. Results are bitwise identical to the old double loop (checked for 1–7 events, 1–3 cycles, n = 1 and 3, with locked patches). Tests: `TimeIntegrationTest/cyclesChainFromTheEndState`, `everyEventIsAppliedOnce`, `badEventTimesAreRejected`.*~~
- [x] ~~**The megathrust mesh has one patch fewer than asked for, and the deepest patch is twice as wide.** `flt2flt` merges the last two patches when the remaining length is exactly `2*W0`, because the floating-point `Wt > Wc/2` comparison is true for equal values. So `megathrust2d.seg` (W0 = W/50) gives 49 patches, and `create_megathrust(..., Nfault)` gives `Nfault-1`. → Use a tolerance in the comparison, or compute the patch count first. Rebuild `kernels/*.mat` afterwards (§6).~~ *Won't fix (2026-11-20): a round-off edge case with no effect beyond one coarser patch at the deepest 8 km. Everything downstream uses `rcv.N`/`rcv.W`, not `Nfault`, and the cached kernels are built on this 49-patch mesh. Fixing it would change the example mesh and force a kernel rebuild. Revisit (one-line tolerance in `flt2flt`) whenever `kernels/*.mat` are rebuilt anyway. Test: `GeometryTest/uniformSegmentMergesTheLastTwoPatches` pins the current mesh (no longer `KnownIssue`).*
- [x] ~~**`geometry.computeFaultDisplacementKernels` nudges a point that sits on a patch centre along the wrong normal.** It uses `rcv.nv(r==0,:)`, where `r==0` is a mask over observation points, not patches. The point gets another patch's normal, or the call errors ("true value outside of the array bounds") when the point's row number is greater than `rcv.N`. → Use `rcv.nv(k,:)`. Also, `delta = 1e-9*mean(ox(:,1))` is 0 or negative when the mean x is ≤ 0. *Fixed: a point on the centre of patch k now moves by `1e-9*rcv.W(k)` along `rcv.nv(k,:)`, so it always gets the value on the +nv side. No current result changes. The only callers that put points on patch centres are the BEM self-problems (`boundary` on `boundary.xc`), where the old indexing happened to pick the right normal and the step was positive. Every BEM kernel of the example (`evl`, `devl`, the hinge kernels) matches the old values to 1e-8 or better, so `kernels/*.mat` don't need regenerating. Tests: `GeometryTest/faultDisplacementKernelOnPatchCentre` (now also checks the value) and `faultDisplacementKernelOnPatchCentreNegativeX` (new). Both fail on the old code.*~~
- [x] ~~**The BEM stress kernels depend on the boundary velocity.** `computeFaultTractionKernelsBem` and `computeShzStressKernelsBem` put `boundary.Vx`/`Vz` into the Dirichlet data of every unit-source column. The displacement-kernel routines use `0.*Vx`. With `bcx`/`bcz ≠ 0` in the boundary `.seg` file (for example `boundarytight2d.seg`), every kernel column picks up a constant term, so the kernels are no longer linear in the source. This has no effect on the current example, because `boundary2d.seg` has zero boundary velocity. → Use `0.*Vx` there too, and treat a prescribed boundary motion as a separate loading term. Test: `KernelsTest/stressKernelsIgnoreBoundaryVelocity`. *Fixed: both routines now use zero Dirichlet data, so all BEM kernels (stress, traction and displacement) solve the homogeneous problem and ignore `boundary.Vx`/`Vz`. The `bcx`/`bcz` columns stay in the `.seg` format for later. The background rate comes from the long-term solution. No current result changes: `boundary2d.seg` has zero boundary velocity, and the example `KK` matches `kernels/evl_orig.mat` to 3e-14. Test: `KernelsTest/stressKernelsIgnoreBoundaryVelocity` (no longer `KnownIssue`; now also checks the mantle kernels).*~~
- [ ] **Future: particular solution for a prescribed boundary motion.** To use `bcx`/`bcz ≠ 0` (for example `boundarytight2d.seg`), add one BEM solve with the prescribed `Vx`, `Vz` and no source, and superpose its stresses and surface velocities once on the cycle response (a loading term in the ODE, not part of the kernels). Not needed for spun-up cycle perturbations.
- [x] ~~**The n = 3 corner flow runs `fsolve` on every call, and prints its output each time.** It also never checked convergence: for dips below about 8° the wedge solve failed silently and returned velocities off by up to 0.6 of the plate rate. *Fixed: the four coefficients now come from one scalar equation for the phase (`fzero`, every root bracketed), with the other three solved linearly; the result is cached per dip and an error is raised if no root meets the boundary conditions. No Optimization Toolbox needed. Tests: `LongTermFlowTest/powerLawCornerFlowBoundaryConditions` and `cornerFlowSolvesStokes` over dips 3–80°, `stokesCheckRejectsWrongRheology`, `powerLawCornerFlowIsQuiet`, `powerLawCoefficientsAreSolvedOncePerDip`.*~~
- [x] ~~**`computeShzStressKernelsBem` added the box correction to the wrong stress components.** For mantle receivers, the half-space stress is stacked as (σxx, σxz, σzz), but the BEM box correction from `geometry.computeFaultStressKernels` was stacked as (σxx, σzz, σxz). So the σxz and σzz corrections were swapped in `LL`. The full-kernel reciprocity test missed it because the correction is only \~1% of `LL`; the correction alone was 65% asymmetric under Betti's theorem. *Fixed in `2eede92` by reordering the rows, and `kernels/evl_orig.mat` was regenerated. `LL` changes by 0.5%, and `KK`, `KL` and `LK` are unchanged. The growing modes and the fastest decay rates don't move, but the slow modes do: in `MAIN_analyticsolve_singleeqcycle.m`, strain rates 50 years after the earthquake change by up to \~20%, and late-interseismic surface velocities by up to 0.04 Vpl. Results computed before the fix carry this error. Test: `BoundaryCorrectionTest/boxCorrectionIsReciprocal` (correction alone, 0.15% asymmetric after the fix).*~~

### Affects unused or legacy code (delete rather than fix, see §1-§3)

- [x] ~~**`computeStressPlaneStrainTriangleShearZone.m` (analytic, 2,238 lines) gives wrong stresses.** It disagrees with the finite-difference version everywhere, including the traction-free surface. Line 84 overwrites the summed `u2d2` with a single term, which looks like leftover debugging. The eigenstrain is then removed using `e22`/`e23`/`e33` after they have been overwritten with the total strain. And `Omega(x2,x3)*e22` is a matrix product, so vector input errors. → Delete it. The finite-difference version is correct: it passes the surface-traction and reciprocity tests. Tests: `GreensFunctionsTest/triangleStressAnalytic*`. *Deleted, with its two tests and its docs entry. Nothing ever called it.*~~
- [x] ~~`computeStressFault_bem.m` errors with "Undefined function `computeFullTractionKernels`". *Deleted, with its docs entry and its test.* Test: `LegacyCodeTest/stressFaultBemRuns`.~~
- [x] ~~`getStrainratesLongterm_old.m` declares `getStrainratesLongterm`. *Deleted.* Test: `LegacyCodeTest/functionNamesMatchFileNames`.~~
- [x] ~~Four file names exist in both `functions/` and `+geometry/`. *Fixed: the `functions/` copies are deleted.* Test: `LegacyCodeTest/noDuplicateFileNames`.~~
- [x] ~~`reconditionODEoperator.m`: the help text says growing eigenvalues are replaced by `-eps`, but the code uses `-0.1` s⁻¹, a 10-second relaxation. `computePeriodicInitialCondition.m` prints the condition number on every call. Tests: `SolversTest/reconditionReplacesGrowingModesByMinusEps`, `SolversTest/legacyPeriodicSolverIsQuiet`.~~ *Deleted (2026-11-21, §3).*
- [x] ~~`pcolor_logvals.m` draws nothing. *Deleted, with its docs entry and its test.* Test: `UtilitiesTest/pcolorLogValsDrawsSomething`.~~
- [x] ~~`create_megathrust.m` writes `megathrust2d.seg` into the current folder. Test: `UtilitiesTest/createMegathrustLeavesWorkingFolderClean`.~~ *Fixed (2026-11-21): renamed `createMegathrust`, writes to a temporary file.*

### Discretisation: growing modes of the mantle kernel

- [x] ~~The collocated mantle-mantle block `LL` has eigenvalues with positive real part: 43 on the example mesh (max 2.4e3), and some even on a 24-triangle mesh. They appear on every mesh tried, including structured uniform meshes, so remeshing is not a fix.~~ *Resolved by §0a: the default kernels are cell-averaged and certified and have none. `negateGrowingModes` (semi-analytical solvers) and `kernelmodify` (numerical solves) are no longer the default treatment: the solvers no longer call `negateGrowingModes` and raise an error on growing modes; both remain as opt-in fixes for centroid kernels (`'cellaverage', 1`).* Test: `KernelsTest/cachedExampleKernelIsCertified` (was `unmodifiedMantleKernelHasNoGrowingModes`).
- [x] ~~Added a damped kernel for numerical solves: `kernelmodify = 2` / `modifyGrowingModes(evl, 2)` sets the growing eigenvalues to -max|lambda| instead of 0. Zeroed modes are neutral and can accumulate over many cycles (reported to run away within 30 cycles at 1e18 Pa s for n = 1, and to need ~100 cycles to settle for n = 3). Tests: `KernelsTest/dampedKernelHasOnlyDecayingModes`, `modifiedKernelsHaveNoGrowingModesForPowerLaw`.~~ *Kept as an option for centroid kernels; not needed with the default certified kernels (§0a).*
- [x] ~~Decide the default for numerical solves.~~ *Decided: the default kernels are cell-averaged and certified (§0a); `kernelmodify` now defaults to 0 and is kept as the earlier empirical fix for centroid kernels.*

### Confirmed correct (safe to rely on during the clean-up)

- The kernel sign conventions agree with each other: KK, KL/LK and LL satisfy Betti reciprocity (ratio ≈ +1) on both the example and the synthetic mesh. This full-kernel check is dominated by the half-space part. The BEM box correction is now checked separately by `BoundaryCorrectionTest`, which caught the swapped rows above. `computeFaultStressKernels`' `-Sxz` convention matches `computeFaultTractionKernels`.
- **Old corner flow = new corner flow** (n = 1) to 1e-14. So `compute_cornerflow_velocityfield` could be replaced by `computeCornerflowVelocityfield(…, 1)` in `computeLongTermVelocity` without changing results (done in §1; the old routine is deleted). Both satisfy their boundary conditions and are divergence-free.
- **Duplicate copies give identical results**: `functions/LDdispHS`, `EdgeStress`, `computedisplacementGFs`, `computetractionkernels`, and the old edge-dislocation traction kernels, which match the line-dislocation kernels. The deletions in §2 change no numbers.
- `computeSiteVelocitiesBem` gives the same result as `devl` from `computeAllDisplacementKernelsBem`, to 1e-12.
- **The ODE path (`ode_viscoelastic` + `ode45`) matches `analyticSolveIVP`** for a locked fault with a Newtonian mantle, through a coseismic step. The solvers match `expm` to round-off.
- **The spun-up numerical cycle matches the semi-analytical periodic cycle** (n = 1, fully locked fault, same operator: unmodified `LL` with growing modes negated). `runImposedEarthquakeCycles` converges to `computePeriodicInitialConditionBalanced` + `analyticSolveIVP` at the rate of the slowest mode (gap shrinks ~250x per two cycles for rho = 0.063) and agrees to 1e-10 after 8 cycles; one cycle alone is still 2e-3 off. Test: `TimeIntegrationTest/spunUpLockedNewtonianCycleMatchesAnalyticPeriodicCycle`. Exact only for a locked fault: the numerical fault is rate-strengthening (`d log v/dt = (dtau/dt)/Asigma`), the analytic one linear viscous. The same comparison on the default certified kernel, used unmodified in both paths (no growing modes, so `negateGrowingModes` is a no-op), agrees to 1.3e-10 after 15 cycles (rho = 0.26): `TimeIntegrationTest/spunUpCycleOnCertifiedKernelMatchesAnalyticPeriodicCycle`.

## 0a. Galerkin (cell-averaged), certified mantle kernels: no growing modes at the source (adopted as the default, 2026-11-20)

**Goal.** Build a mantle-mantle kernel that has no growing modes, so that neither `kernelmodify` nor `negateGrowingModes` is needed, for linear (Maxwell) and power-law (n = 3) rheology alike, in both the numerical and the semi-analytical solvers. The example mesh in `inputs/` stays as it is: only the kernels change.

**Diagnosis.**
- The continuous operator from a deviatoric eigenstrain to the stress it causes is symmetric and negative semidefinite: the elastic energy of any eigenstrain distribution is non-negative, so stress always opposes the eigenstrain. A correct discretisation keeps both properties and cannot have growing modes.
- I2EC uses piecewise-constant strain per triangle (the source integral is exact, Barbot 2018), but samples the receiver stress at the **centroid** (collocation, i.e. a point test function). For neighbouring cells the stress varies strongly across the receiving cell, so the point value is not the energy-consistent driving stress. This breaks the symmetry (area-weighted `LL` is 15.6% asymmetric on the example mesh, 13-14% for neighbours, 0.5% beyond 4 cell sizes) and the sign: the symmetric part has 65 positive eigenvalues.
- Using the **cell average** of the receiver stress instead (Galerkin with P0 elements: same piecewise-constant function for trial strain and test stress) keeps symmetry and sign. The average is computed with a 3-point (degree 2) or 7-point (Dunavant degree 5) interior quadrature rule, so it is approximately Galerkin (adjacent pairs have edge discontinuities and corner log-singularities a fixed rule cannot integrate exactly).
- Ruled out: the boundary correction (the half-space kernel alone has the same 43 modes) and the pressure removal in the deviatoric projection (identical in the centroid and averaged kernels; for a deviatoric eigenstrain only the deviatoric stress does work, so dropping the pressure cannot make the operator indefinite). Collocation asymmetry alone is not the cause either: symmetrising the centroid kernel gives 65 growing modes, not 0.

**Results so far.**

| Kernel (example mesh, 521 cells) | Growing modes | max Re(lambda)/max abs(lambda) |
|---|---|---|
| centroid, half-space only | 43 | 5.7e-2 |
| centroid, full BEM (cached `evl_orig`) | 43 | 5.7e-2 |
| cell average 3-point, half-space only | 0 | -2.2e-3 |
| cell average 7-point, half-space only | 0 | -2.2e-3 |
| cell average 3-point, full BEM | 0 | -2.2e-3 |

- Averaging changes `LL` entries by a median 17-31% for cells within 2 cell sizes, 7% at 2-4, 1% at 4-10 and 0.2% beyond: only the near field changes.
- First toy mesh (PDE Toolbox, Hmax 40 km, 127 cells): 8 growing modes with centroid, 0 with 3- and 7-point averaging. Structured diagonal-split mesh, 0.8 km cells, 1250 cells (half-space only, circular-inclusion setup): 244 growing modes with centroid, 0 with 3-point averaging. 24-triangle test mesh: growing modes with centroid, none with averaging (unit test).
- Cost: about 3x (3-point) the time of the `KL` and `LL` blocks; 284 s for the full example kernel.

**Toy meshes (full BEM, example boundary).** Growing modes, centroid → 3-point → 7-point:

| Mesh | Cells | Centroid | 3-point | 7-point | max Re/max abs (3-point) |
|---|---|---|---|---|---|
| PDE Toolbox outline, Hmax 40 km | 127 | 8 | 0 | 0 | -9.7e-3 |
| PDE Toolbox outline, Hmax 25 km | 320 | 25 | 0 | 0 | -3.4e-3 |
| PDE Toolbox outline, Hmax 15 km | 876 | 78 | 0 | – | -1.1e-3 |
| structured diagonal, 25 km, box z -150 to -50 km | 96 | 41 | 29 | 35 | +6.5e-1 |
| structured diagonal, 15 km, same box | 280 | 114 | 70 | 70 | +9.5e-1 |
| structured diagonal, 10 km, same box | 600 | 208 | 134 | – | +5.4e-1 |
| structured criss-cross, 25 / 15 / 10 km, same box | 192 / 560 / 1200 | 25 / 97 / 223 | 12 / 42 / 113 | 7 / – / – | ~1e-7 |
| structured diagonal, 2.5 km, box z -120 to -100 km | 128 | 18 | 0 | – | -7.8e-3 |
| structured diagonal, 0.8 km, box z -120 to -100 km (half-space only) | 1250 | 244 | 0 | – | – |
| 24-triangle test mesh | 24 | 2 | 0 | 0 | – |

- Unstructured meshes: averaging removes every growing mode.
- **Shallow structured meshes still fail, and differently.** Even the centroid kernel has a growing eigenvalue as large as the largest eigenvalue (max Re/max abs up to 1.0), which is not the subtle collocation effect. The BEM boundary is not the cause (half-space-only kernel matches the BEM kernel to 1e-9, same counts). The triangle-stress function reports hundreds of "points corrected with double-exponential integration" on these meshes. Hypothesis under test: when the grid is aligned with the surface and starts at a whole multiple of the cell size, centroids and quadrature points lie exactly on the extensions of the *image* triangles' edges, where the closed-form half-space solution is singular; only non-finite values are corrected, so finite-but-wrong values (e.g. an angle branch jump) pass through. The image is far for deep boxes (which pass) and close for shallow ones (which fail). **Confirmed by shifting the same 15 km mesh by a fraction of a cell** (half-space kernel):

  | Shallow 15 km structured mesh | Centroid | 3-point average |
  |---|---|---|
  | aligned (as in the table) | 112 growing, max Re/max abs 1.0, area-weighted asymmetry 141% | 72 growing, 0.94 |
  | shifted 0.37 cell | 46 growing, 0.14, asymmetry 15% | **0**, asymmetry 1.5% |
  | shifted 0.5 cell | 46 growing, 0.14 | **0** |

  Once nothing is aligned, the mesh behaves like all the others and averaging removes every growing mode. The aligned mesh has a few grossly wrong kernel entries (141% asymmetry) that are finite, so the function's non-finite check does not catch them. This is a separate bug in `computeStressPlaneStrainTriangleShearZoneFiniteDifference` and may also explain growing modes seen on structured meshes in earlier tests. **Found and fixed.** The bad entries are ~9e10 (typical off-diagonal ~90, self term ~-1.6e4). Every one has its receiver centroid exactly on the line through an edge of the *image* triangle (distance 0). The closed-form displacement contains `atan(N/D)` with `D = (x3+m3)*a2 - a3*(m2-x2)`, which is zero on that line; `atan` jumps by pi there although the field is smooth, and the 1e-6 finite-difference stencil across the jump gives a finite but enormous strain. The same denominator appears for lines through the source's own edges.
- [x] ~~Fix the triangle-stress evaluation at aligned points.~~ *Fixed in `computeStressPlaneStrainTriangleShearZoneFiniteDifference` (d5279e7): points within 1e-5 l of a line through a source or image edge (l = shortest edge; not on the triangle's own edges or near its vertices) are evaluated as the mean of two points shifted +-1e-3 l normal to the line (one-sided at the surface); other points are unchanged. On the aligned shallow 15 km mesh the centroid kernel now has 46 growing modes (max Re/max abs 0.14, as for the shifted mesh) and the 3-point kernel 0 (asymmetry 1.5%). Test: `GreensFunctionsTest/triangleStressIsSmoothOnEdgeLines` (points on an image-edge line and on a source-edge line outside the triangle; fails on the old code).*
- [x] ~~Rerun the structured toy meshes with the fix.~~ Growing modes, centroid → 3-point: diagonal 25 / 15 / 10 km: 12 → 0, 46 → 0, 108 → 0 (max Re/max abs -8e-3, -3e-3, -1e-3). Criss-cross 25 / 15 / 10 km: 25 → 12, 97 → 42, 223 → 113, all at ~1e-7 of the spectrum: zero-energy modes (compatible, divergence-free piecewise-constant strains that criss-cross meshes admit) perturbed by quadrature. Removed by the certification below.
- [x] ~~Check whether the example mesh or the cached kernels have any points on such lines.~~ *No effect on the example: `KernelsTest/cachedStressKernelsAreReproducible` passes against the regenerated `kernels/evl_orig.mat`, which certifies with nothing clipped.*
- Criss-cross meshes: the remaining modes are at ~1e-7 of the spectrum. Criss-cross meshes admit compatible, divergence-free piecewise-constant strains that produce no stress (true zero-energy modes); the approximate quadrature perturbs these zeros slightly. *Resolved: the certification sets them to exactly zero (see below).*

**Linear comparison (example mesh, periodic cycle, 500-yr recurrence).** Centroid kernel + `negateGrowingModes` vs 3-point cell-averaged kernel with no modification:
- The cell-averaged operator `A` has no growing modes; `negateGrowingModes` leaves it unchanged (relative change 1e-13). With the centroid kernel it has 43 and negation changes `A` by 3.6e-5 (viscosity 1e19 / 10^18.5 Pa s) or 3.6e-4 (10x lower).
- Oceanic 1e19, wedge 10^18.5 Pa s: differences (relative, whole field) through the cycle are ≤ 1.7% in fault slip rate, ≤ 1.9% in surface vx, ≤ 3.2% in surface vz, and 3-15% in mantle strain rate.
- Oceanic 1e18, wedge 10^17.5 Pa s: up to 26% in fault slip rate (at 20 yr), 15% in surface vx, 36% in surface vz (late in the cycle) and 46-88% in mantle strain rate (5-20 yr).
- These differences combine the near-field kernel change and the effect of negation; neither run is a converged reference. To tell which is closer to the truth, compare both against a refined mesh (or the circular-inclusion test for the kernel itself).

**Implemented (d5279e7).**
- `computeAllStressKernelsBem(..., 'cellaverage', 1 | 3 | 7)`: 1 = centroid (the earlier kernel), 3 (now the default) or 7 = cell average for the two blocks with a mantle receiver (`KL`, `LL`). `KK` and `LK` are unchanged (fault receivers). Averaging `KL` also makes it consistent with `LK`, which already treats each mantle cell as a whole source.
- Implementation: the BEM receiver code reads only `shz.xc`, so the kernel is evaluated at each quadrature point in turn (temporarily moving `shz.xc`, restored by `onCleanup`) and the weighted sum is taken.
- Tests (`KernelsTest`): `cellAveragedMantleKernelHasNoGrowingModes` (3 and 7 points, also checks area-weighted symmetry < 2% and that centroids are restored), `centroidCollocationHasGrowingModesOnTestMesh`, `defaultIsCertifiedThreePointAverage`, invalid values rejected.

**Checks before adopting it (all done).**
- [x] ~~Toy meshes (used only to test that averaging removes the growing modes on any mesh, not to replace the example mesh).~~ See the tables above: passes on every mesh once the edge-line bug was fixed; criss-cross meshes after certification.
- [x] ~~Linear comparison: periodic cycle on the example mesh, centroid kernel + `negateGrowingModes` vs cell-averaged kernel with no modification, two viscosity levels. Report differences in fault slip rate, mantle strain rate and surface velocity through the cycle; check that the cell-averaged operator `A` has no growing modes (so negation does nothing).~~ Done, see above; differences need a converged reference to interpret.
- [x] ~~Circular-inclusion accuracy test: structured 0.8 km mesh at 110 km depth, sharp and tapered sources, diagonal transect, against the exact circular inclusion (Eshelby; interior value -mu*eps/(2(1-nu))). Compares centroid, `kernelmodify` 1 and 2, and 3-point averaging; repeats the test behind `Figure_regkernels_comparison.pdf` (which used a separate kernel library) with I2EC's triangle kernels.~~ Done. Growing modes: centroid 244, zeroed 0, damped 0, 3-point average 0. The reference matches Eshelby (interior -0.6658 vs -0.6667 mu*eps). RMS error of the sharp-source stress along the transect (units of mu*eps; collocated kernels scored against exact point values, the averaged kernel against exact cell averages):

  | Kernel | inside (r < 3 km) | edge (3-5.5 km) | outside (r > 5.5 km) |
  |---|---|---|---|
  | centroid (original) | 0.003 | 0.018 | 0.002 |
  | `kernelmodify = 1` (zeroed) | 0.004 | 0.018 | 0.003 |
  | `kernelmodify = 2` (damped) | 0.047 | 0.050 | 0.024 |
  | cell average, 3-point | 0.004 | 0.016 | 0.002 |

  The cell-averaged kernel is as accurate as the original (for the quantity it represents) and has no growing modes. Zeroing changes this smooth-ish source by < 1%; damping is 10-20x less accurate. On this mesh (source radius = 5 cells) all kernels agree to a few % along the transect, unlike the much larger changes in `Figure_regkernels_comparison.pdf`, which used a different kernel library and a finer square-cell mesh.
- [x] ~~Power-law check: linearise the n = 3 mantle system about random stress states with the cell-averaged kernel and confirm no growing modes.~~ Done on the example mesh (3-point full BEM kernel): 0 growing modes for n = 1 and n = 3 in every trial; all eigenvalues strictly negative (max Re/max abs -1e-3 for n = 1, -6e-7 for n = 3).
- [x] ~~Rebuild the example kernels with averaging (same mesh in `inputs/`; only the kernels are recomputed) and compare `MAIN_analyticsolve_singleeqcycle.m` and `MAIN_i2ec_singleeqcycle.m` outputs with the current ones.~~ Done (2026-11-20, scratch copies of both drivers; the repo scripts and `kernels/` unchanged).
  - Certified 3-point example kernel: built in 165 s; certificate: asymmetry removed 1.35%, largest eigenvalue -8.3e-4 of the spectrum, **0 eigenvalues clipped** (on the example mesh the certification only removes quadrature asymmetry).
  - Analytic driver (centroid kernel + `negateGrowingModes` vs certified kernel, negation then a no-op): state 4e-5 to 1.8% through the cycle; late-interseismic surface vx 1.0%, vz 2.6%; surface displacement time series ux ≤ 1.2%, uz mostly ≤ 6% (up to 14% at a few points where uz is near zero).
  - Numerical driver (centroid kernel + `kernelmodify = 1` vs certified kernel, no modification), last of 5 cycles: surface vx median 0.1%, max 0.8%; vz median 1.0%, max 2.0%; fault slip rate max 0.7%; mantle strain rate median 1.9%, max 4.2%.
  - **The certified kernel ran 3.6x faster in `ode45`**: 239 s vs 850 s for the 5 cycles (fewer, larger steps).
  - Figure: surface velocity profiles from both drivers overlap at plot scale.

**Certification (d5279e7; decided 2026-11-20).**
- Decision: use 3-point averaging, for the mantle-mantle block only. The fault blocks stay collocated: for the straight, non-interacting faults in this repo the fault kernels are negative semidefinite already, so Galerkin is not needed there.
- `functions/certifyMantleKernel.m`, called by `computeAllStressKernelsBem` whenever `'cellaverage'` > 1: symmetrises the area-weighted `LL` (removes only quadrature error), checks that the largest eigenvalue of the symmetric part is ≤ tolerance × its largest magnitude (error otherwise), sets eigenvalues in (0, tolerance] to exactly 0, and stores `evl.certificate` (asymmetry removed, max eigenvalue, number clipped, tolerance).
- Guarantee: a symmetric negative semidefinite area-weighted kernel gives no growing modes for any viscosity distribution and any power-law exponent at any stress state (the linearised system L*D is similar to X*S*X with X symmetric positive definite).
- Tolerance 1e-3 (relative to the largest eigenvalue magnitude), from the measured margin: quadrature-perturbed zero-energy modes reach 1.4e-4 (criss-cross 25 km, 3-point); the centroid example kernel, which is truly indefinite, has 65 positive eigenvalues up to 4.3e-2; the 3-point example kernel has none (largest -8.3e-4). The quadrature error itself is ~1e-2 (antisymmetric part of the example kernel 8.9e-3 of the spectrum), so clipping at 1e-3 is below it.
- Tests (`KernelsTest`): certified kernel exactly symmetric and negative semidefinite (3 and 7 points), `certificationRejectsCentroidKernel` (example centroid kernel errors), `certifiedKernelHasNoGrowingModesForPowerLaw` (n = 1 and 3, random stress states over 4 orders of magnitude and prefactors over 4 orders).
- [x] ~~Confirm the certification passes on the criss-cross meshes at 1e-3.~~ Criss-cross 25 km: asymmetry 2.3%, largest symmetric eigenvalue 1.4e-4, 20 clipped; 15 km: 2.3%, 1.6e-4, 90 clipped. Certified kernels: largest real part ~1e-15 of the spectrum (round-off), all eigenvalues real. (10 km, 1200 cells, not rerun; same margin expected.)

**Bug found by the certification: nearly axis-aligned edges (2026-11-20).** The default kernel of a finer PDE Toolbox mesh (`createShearZoneMesh('Hmax', 10)`, 1955 cells) failed certification with a positive eigenvalue at 0.997 of the spectrum. 17 cells had self terms up to 130x too large, some positive: equilateral cells in regular rows whose "horizontal" base differed in depth by ~1e-9 m (floating-point noise from the mesher; edge direction [1, -9.9e-14]). The closed-form displacement (`computeDisplacementPlaneStrainTriangleShearZone`) is exact for an exactly horizontal or vertical edge but divides by the edge's small direction component when it is only nearly aligned: u2 error per unit strain 66 at a 1e-10 m tilt of a 10 km edge, 0.15 at 1e-8 m, 8e-6 at 1e-3 m (6e-10 when exactly horizontal).
- [x] ~~Fixed (4ed5e8d): edges horizontal or vertical to within 1e-7 of their length are snapped to exactly horizontal or vertical before evaluation (vertices move by at most 1e-7 of an edge). All 1955 self terms of that mesh are now normal (-2.2e4 to -1.8e4, median -2.0e4). Test: `GreensFunctionsTest/triangleDisplacementWithNearlyAlignedEdges` (fails on the old code).~~
- [x] ~~Confirm the full Hmax 10 kernel now certifies.~~ It does: built in 765 s; asymmetry removed 0.89%, largest eigenvalue -3.8e-4 of the spectrum, nothing clipped.

**Late-interseismic strain rates at the Maxwell-time endmembers (n = 1, periodic = fully spun-up cycle, example mesh, uniform viscosity, t = 0.98 T).** Centroid kernel + `kernelmodify = 1` vs certified 3-point kernel:

| Maxwell time / T_recur | Kernel | Operator modes with Re ≥ 0 | Median strain rate / long-term |
|---|---|---|---|
| 1000 | `kernelmodify = 1` | 25 slightly positive, 42 exactly zero | 1.09 |
| 1000 | certified | none | 1.00 |
| 0.001 | `kernelmodify = 1` | 27 slightly positive, 30 exactly zero | 11.7 |
| 0.001 | certified | none (slowest mode e-folds in 0.35 T) | 4.3 |

The zeroed modes are neutral, so in the spun-up cycle they hold cell-scale noise (speckled strain-rate maps). With a weak mantle (Maxwell time 0.001 T) the noise dominates: field difference 42%, per-cell ratio certified/zeroed 0.09-1.9 (5th-95th percentile), fault slip rate 5%. With a strong mantle (1000 T) the difference is confined to the slowly deforming regions (field difference 2%, per-cell ratio 0.17-1.4). The certified fields are smooth.

**Adopted (d5279e7; solvers stop negating in 7d0e61f).**
- [x] ~~`computeAllStressKernelsBem` defaults: `'cellaverage', 3`, `'kernelcertify', true` (new option), `'kernelmodify', 0`.~~ `kernelcertify` with `'cellaverage', 1` is an error (centroid kernels cannot be certified); the earlier centroid kernel is `'cellaverage', 1, 'kernelcertify', false` (plus `'kernelmodify', 1` for the old numerical default).
- [x] ~~`kernelmodify` and `modifyGrowingModes` kept as the earlier, empirical fix for centroid kernels (not removed).~~
- [x] ~~`kernels/evl_orig.mat` regenerated with the default (certified) kernel; it now has an `evl_orig.certificate` field (asymmetry 1.35%, max eigenvalue -8.3e-4, 0 clipped). `devl.mat` is unchanged (displacement kernels do not depend on the mantle receivers).~~
- [x] ~~Help text updated: `computeAllStressKernelsBem`, `modifyGrowingModes`, `ode_viscoelastic`, `runImposedEarthquakeCycles`, `assembleStressKernel`, `computeInterseismicSteadystate`; comment in `MAIN_analyticsolve_singleeqcycle.m`. `docs/theory.md` growing-mode section rewritten (prevention by default; earlier fixes for centroid kernels); `docs/inputs.md` describes the cached kernel.~~
- [x] ~~Tests: `KernelsTest` uses an explicit centroid kernel where the test is about centroid kernels, plus `defaultIsCertifiedThreePointAverage` and `cachedExampleKernelIsCertified` (replaces the `KnownIssue` test `unmodifiedMantleKernelHasNoGrowingModes`); `TimeIntegrationTest/spunUpLockedNewtonianCycleMatchesAnalyticPeriodicCycle` pinned to the centroid kernel it was written for (its last check needs a slow mode).~~
- [x] ~~Negation is no longer called once certification passes (2026-11-20).~~ `analyticSolveIVP` and `computePeriodicInitialConditionBalanced` use a new `eigBalanced` (balanced eigendecomposition, shared with `negateGrowingModes`) and raise `I2EC:growingModes` if the operator has a mode with Re(lambda) > 1e-10 max|lambda|, instead of negating it. `MAIN_analyticsolve_singleeqcycle.m` no longer calls `negateGrowingModes`. `negateGrowingModes` stays as an opt-in for centroid kernels: negate the operator explicitly and build the forcing with it. Tests: `SolversTest/solversRejectGrowingModes`, `eigBalancedMatchesEig`; `test_negate_exact_solver.m` test 7 passes the explicitly negated operator.
- [x] ~~Examples and docs highlight the default: both MAIN scripts explain the default kernels (3-point cell average, certified) and show the call for the earlier centroid kernels with `kernelmodify`; `docs/quickstart.md` has a "Kernel options" table (default, 7-point, centroid, centroid + `kernelmodify` 1 or 2, explicit negation); `docs/theory.md` marks negation as opt-in and the treatment caveats as centroid-only, and its practical guidance points to `evl.certificate`; `README.md` lists certified mantle kernels as a feature.~~
- [ ] Optionally rename `evl_orig` (it is no longer an "original, unmodified" kernel in the old sense).
- [x] ~~**Certify the coupled fault-mantle operator, not only `LL`.**~~ `certifyMantleKernel` guarantees no growing modes for the mantle alone (fault locked or held fixed), but the solvers use the full operator `[KK LK; KL LL]` with a creeping fault, and the fault blocks are not certified (decision above: `KK` alone is negative semidefinite for the straight faults here). The coupling can move eigenvalues even when `KK` and `LL` are each fine, and its two halves are discretised differently (`KL` cell-averaged over each mantle cell, `LK` collocated at fault patch centres). Today this is only checked case by case: the example 2 operator and the linearised operator of the former `MAIN_i2ec_singleeqcycle.m` have no growing modes (2026-11-20). → Weight fault rows by patch width and mantle rows by cell area, `S = diag([W; area; area]) * [KK LK; KL LL]` (symmetric in the continuum by Betti reciprocity), and check that its symmetric part is negative semidefinite, as for `LL`. That guarantees no growing modes for any fault law (linear viscous, or rate-strengthening linearised about any slip rate), any viscosity distribution and any n, with the fault creeping. Report the asymmetry between the area-weighted `KL` and the width-weighted `LK` (it measures how inconsistent the two coupling blocks are), and decide whether to symmetrise the coupling blocks too or only check them. Locked patches drop out of the operator, so certify the full fault (any locked subset is then covered).
  - *Findings so far (2026-11-21, branch `coupled-certification`; scripts in `investigations/coupled_certification/`).*
    - **Weights.** Betti holds with a factor 2 from the deviatoric work pairing (σ22 e22 + σ33 e33 + 2σ23 e23 = 2[(σ22−σ33)/2 · e22 + σ23 e23]): `W .* LK(:,:,c) ≈ 2 * (a .* KL(:,:,c))'` (fitted 1.997). So weight mantle rows by `2*area`: `S = diag([W; 2a; 2a]) * [KK LK; KL LL]`. Check its symmetric part after Jacobi scaling (`D^-1/2 S D^-1/2`, D = |diag S|; congruence keeps the signs) so the fault and mantle blocks are compared on equal terms; unscaled, the `LL` self terms swamp everything.
    - **Example geometry (`coupledExample.m`, `fault_mantle_gap.png`).** There is a gap: the slab mantle is 19.4 km (normal distance) below patches 1-45, i.e. the elastic plate. Patches 46-49 are 17, 13, 9 and 1 km from the mantle; the down-dip tip of patch 49 nearly touches the corner of the mantle wedge (its collocation point is 5 km away; cells there ~7 km, patches 4.1 km).
    - **Example result.** The coupled operator certifies with the default kernels: largest eigenvalue of the scaled symmetric part −4.3e-3, none ≥ 0 (fault block alone −1.3e-2). The KL/LK mismatch per patch (scaled) is 1.5% median but 30% for patch 49; averaging LK over each patch (Gauss points along it) brings it to 0.7% and 8%. `KK`: width-weighted asymmetry 3%, symmetric part negative definite, as expected for this straight fault; no averaging needed there.
    - **Fault at distance d from the mantle (`coupledDistanceSweep.m`: 100 km fault, 2.5 km patches, structured band of 5 km cells parallel to it).** Default blocks (KL 3-pt, LK collocated): mismatch 1.4% / 13% / 34% / 65% (max over patches) at d = 20 / 5 / 2 / 1 km, but the margin is unchanged (−6.0e-3) down to d = 2 km, −3.9e-3 at 1 km, −6.7e-4 at 0.5 km, and the symmetric part turns indefinite at d ≤ 0.1 km (one positive eigenvalue; also with KL 7-pt and LK 5-pt, and 12 at d = 0 with KL 3-pt/LK 5-pt, where K·diag(D) then has 2 growing modes for some random D). So, as expected: the gap keeps the collocated coupling harmless; only cells within about a patch width of the fault matter.
    - **Fix that works at every d, including the fault lying on cell edges.** Compute LK averaged over each patch accurately (Gauss points along the patch; 32 points converge to 3e-9 at d = 0.1 km and 4e-4 at d = 0, 16 points to 1e-13 at d ≥ 1 km) and take KL from it by reciprocity, `KL(:,:,c) = (W'/2 ./ a) .* LK(:,:,c)'`. Then the coupling blocks are exactly consistent and the coupled operator certifies at all d (largest eigenvalue −6.0e-3 at 2 km, −4.8e-3 at 1 km, −1.5e-3 at 0.1 km, −8.4e-4 at d = 0), with `LL` 3-pt as now. Near the fault the computed cell-averaged KL is the inaccurate block (it differs from the Betti KL by 32% at d = 1 km even with 7 points): a dislocation's stress varies on the scale d, which a fixed cell rule cannot resolve, whereas the triangle stress along a patch is smooth enough for Gauss quadrature.
    - **Zero gap at the wedge (`coupledZeroGap.m`).** Moving the wedge's left edge (x = 201 km) onto the fault tip, so the wedge corner touches the down-dip end of patch 49, changes nothing that matters: `LL` certifies as before (max eigenvalue −8.3e-4), the coupled operator's largest scaled eigenvalue is −4.31e-3 (same as with the 1 km gap), and the example 2 cycle solves. Only patch 49's KL/LK mismatch grows (30% → 45%). The gap at the wedge is not needed for stability; what fails is cells lying along the fault over many patches (sweep above, d ≤ 0.1 km).
    - **Decision (2026-11-21): keep standard BEM kernels; no reciprocity-based coupling and no coupled check in the code.** `KK` and `LK` stay collocated at patch centres, `KL` and `LL` cell-averaged (3-point) with `LL` certified. The coupled operator only loses definiteness when mantle cells sit within a small fraction of a cell size of the fault along much of its length (d ≤ 0.1 km here, with 5 km cells); a mesh with d < cell size is poor experiment design, since the kernels are singular functions (already stated in the docs). Guidance: keep mantle cells at least about one cell size from the fault. Optional: add that sentence next to the singular-kernel caveat in the docs. The scripts stay in `investigations/coupled_certification/` as the record.


- [x] ~~Both safeguards raise errors instead of fixing quietly (2026-11-21).~~ *`certifyMantleKernel` no longer clips eigenvalues in (0, tol] to zero: a kernel with any eigenvalue >= 0 fails (`certifyMantleKernel:zeroEnergyModes` below the 1e-3 threshold, `notNegative` above). `computePeriodicInitialConditionBalanced` no longer skips modes with |lambda| < 1e-16 1/s (an absolute cutoff in 1/s, not machine precision). That cutoff was silently wrong for strong mantles: with the example 2 fault and a uniform mantle of Maxwell time 1e4 / 1e5 / 1e6 T it zeroed 16 / 131 / all slow modes, and the periodic mantle state was 58% / 83% / 100% off (against a direct solve, y0 = -(A G) \ (G c + b) with G from an augmented expm). Without any cutoff the modal solve agrees with that direct solve to 1e-8 / 3e-7 / 1e-6 in the mantle. A cutoff relative to max|lambda| would be wrong too (those modes are 5e-12 to 5e-14 of max|lambda|, set by the fault). Now modes with |Re(lambda)|*T <= 1e-12 (relax by less than 1e-12 per cycle; independent of units) raise `I2EC:neutralModes`; at Maxwell time 1e6 T the slowest mode has |lambda|*T = 3e-9. Weak mantles (1e-5 T) are unaffected either way. No example or cached kernel is affected (nothing was ever clipped on them). Tests: `KernelsTest/certificationRejectsZeroEnergyModes`, `SolversTest/periodicSolverRejectsNeutralModes`, `periodicSolverAcceptsWidelySeparatedTimescales`; `test_negate_exact_solver` test 6c now expects the error.*
- [x] ~~The semi-analytical solvers eigendecompose the operator without balancing (2026-11-21).~~ *`eigBalanced` now calls `eig(A, 'nobalance')` (balancing is an option). A = K./eta has eigenvectors conditioned by the geometry (cond(V) ~ 7e5 for any mantle Maxwell time); balancing raised cond(V) to 2e10 for a strong mantle and lost three digits. Periodic mantle state against a direct solve, example 2 fault: 5e-11 to 2e-9 for t_M/T = 1e-5 to 1e9 (balanced: 4e-12 at example 2, but up to 1e-6 at 1e6); `I2EC:neutralModes` at 1e10. A symmetrising scaling was tried and also loses digits for strong mantles. `test_negate_exact_solver` tolerances 1e-10 → 1e-8 (now ~2e-9 to 5e-9). Docs table updated. Test: `SolversTest/eigBalancedOptionalBalancingGivesSameModes`. The name `eigBalanced` is now misleading; see the naming scheme (§4).*

## 1. Broken or likely incorrect code (fix or delete first)

- [x] ~~**`functions/computeStressFault_bem.m` is broken.** *Deleted, with its docs entry and `LegacyCodeTest/stressFaultBemRuns`.* It called `computeFullTractionKernels`, which doesn't exist anywhere in the repo. Nothing called this function, but the docs listed it (`docs/api/kernels.rst`).~~
- [x] ~~**`functions/getStrainratesLongterm_old.m` has the wrong function name inside it.** *Deleted.* The file declares `function ... = getStrainratesLongterm(...)`, so the declared name doesn't match the filename. Nothing calls it. → Delete it.~~
- [x] ~~**`functions/pcolor_logvals.m` is unfinished.** *Deleted, with its docs entry and `UtilitiesTest/pcolorLogValsDrawsSomething`.* It builds `pos_logvals` and `neg_logvals` and then returns without plotting or outputting anything. Nothing calls it, but the docs list it (`docs/api/utilities.rst`). → Finish it or delete it.~~
- [x] ~~**Two corner-flow implementations.** `compute_cornerflow_velocityfield.m` (old, JPL 2023, Newtonian only) is still used by `computeLongTermVelocity.m` (and so by `examples/plotLongtermVelocity.m`). `getStrainratesLongterm.m` uses the newer `computeCornerflowVelocityfield.m` (supports n = 1 and n = 3). The test suite shows the two agree to 1e-14 for n = 1 (`LongTermFlowTest/oldAndNewNewtonianCornerFlowAgree`). → Switch `computeLongTermVelocity` to `computeCornerflowVelocityfield(…, 1)` and delete the old one.~~ *Done: `computeLongTermVelocity` uses `computeCornerflowVelocityfield` (n = 1) and `compute_cornerflow_velocityfield.m` is deleted, with its docs entry and `LongTermFlowTest/oldAndNewNewtonianCornerFlowAgree`. Velocities change by 2e-14 relative on the `plotLongtermVelocity` grid.*
- [x] ~~**`computeLongTermVelocity.m` and `getStrainratesLongterm.m` repeat the same translate/rescale block** (`ocean_translate_params`, `continental_translate_params`, the `- 1` shift). → Pull it into one helper so the two can't drift apart.~~ *Done: new `cornerFlowVelocity(dip, x, z, region, shift, n)` applies the shift (with the 1 m) and picks the region's column; both callers use it. Strain rates change by under 1e-14 relative (n = 1 and 3). Tests: `LongTermFlowTest/regionVelocityIsShiftedCornerFlow`, `unknownRegionIsRejected`.*
- [x] ~~**`computeInterseismicSteadystate.m` hard-codes geometry inputs:** `[0,20e3]` and `[-140e3,35e3]` are passed to `getStrainratesLongterm`. → Make them arguments, or document them. Moot if the function is deleted (§3): since the notebook was removed, nothing in the tracked code calls it.~~ *Moot: deleted.*
- [x] ~~**The deviatoric `LL` kernel's index construction is unverified.** `computeAllStressKernelsBem.m` (around line 200) builds the deviatoric `LL` from the full sigma22/sigma23/sigma33 kernels with e33 = -e22 imposed, and its comment says the indices still need checking. The theory page used to mention this; it now lives here. Indirect evidence that it is right: Betti reciprocity and certification of `LL` pass. → Add a direct test against the deviatoric stress of a single strained cell.~~ *Verified (2026-11-21): the construction is right. Stress of single strained cells, computed directly, matches `evl.LL` and `evl.LK` to 1e-12 (swapping the off-diagonal blocks would be 2% off, which the symmetry checks cannot see). The stale code comment (which also had the index order backwards) is replaced by the correct description. Test: `KernelsTest/deviatoricKernelsMatchDirectStress`.*

## 2. Duplicate code: `+geometry/` vs `functions/`

Several files exist in both places. Package code calls the `geometry.` versions, so the `functions/` copies are leftovers from older versions.

| File | `+geometry/` copy | `functions/` copy | Action |
|---|---|---|---|
| `LDdispHS.m` | used (`geometry.LDdispHS`) | ~~identical, unused~~ **deleted** | done |
| `EdgeStress.m` | ~~only used by the unused `geometry.computetractionkernels`~~ **deleted** | ~~identical, unused~~ **deleted** | done |
| `computedisplacementGFs.m` | ~~unused~~ **deleted** | ~~unused~~ **deleted** | done (replaced by `computeFaultDisplacementKernels`) |
| `computetractionkernels.m` | ~~unused~~ **deleted** | ~~unused~~ **deleted** | done (replaced by `computeFaultTractionKernels`) |

- [x] ~~Delete the `functions/` copies of the four files above.~~
- [x] ~~Delete the unused `+geometry/` versions: `computedisplacementGFs`, `computetractionkernels` and `EdgeStress` (only `computetractionkernels` used it). *Deleted, with their equivalence tests in `LegacyCodeTest`. Checked first: no caller in the repo, the local `expt*`/`experiments` folders or the git history since 2023, and the current kernels reproduce them to round-off for the example's megathrust, box and hinge, the small test model, and 14 synthetic faults (dips 5–90°, surface-breaking and buried). The current kernels also cover opening, which these did not.*~~
- [x] ~~Remove the docs mentions of the `functions/` copies (the note in `docs/api/kernels.rst`; they had no API entries of their own).~~
- [x] ~~Decide on one rule for what goes where. For example: low-level Green's functions and geometry classes go in `+geometry/`, and the `*Bem` assembly, solvers and plotting go in `functions/`. Write the rule down in the README or docs.~~ *Done: in the naming rules (§4).*

## 3. Unused functions in `functions/`

No caller in the repo. Each one needs a keep-or-delete decision.

- [x] ~~`computeKLmodesMesh.m`: KL/eigenmode basis, looks experimental. Also, the help text calls it `compute_eigmodes_mesh` (out-of-date name).~~ *Deleted (2026-11-21).*
- [x] ~~`computePeriodicInitialCondition.m` (with `reconditionODEoperator.m`, which only it uses): replaced by `computePeriodicInitialConditionBalanced.m` and `negateGrowingModes.m`. → Delete both, or keep one as a documented fallback.~~ *Both deleted (2026-11-21); `computePeriodicInitialConditionBalanced` was renamed `computePeriodicInitialCondition`.*
- [x] ~~`computeSiteVelocitiesBem.m`: nothing calls it, and it carries a block of commented-out `devl`-based code.~~ *Deleted (2026-11-21).*
- [x] ~~`computeStressPlaneStrainTriangleShearZone.m` (**2,238 lines**, Barbot 2018 analytic stress): the code actually uses `computeStressPlaneStrainTriangleShearZoneFiniteDifference.m` instead. The tests show it gives **wrong stresses** (§0), so it cannot serve as a reference. → Delete it. This is the largest single file in the repo. *Deleted.*~~
- [x] ~~`computeInterseismicSteadystate.m`: its only caller was the notebook, now deleted (the remaining callers are in the git-ignored `experiments/`). Interseismic steady state with a pinned fault and a Laplacian-smoothed mantle; the semi-analytical driver does not use it. → Delete, or keep and add it to a driver script.~~ *Deleted, with `SolversTest/interseismicSteadyStateRuns` and its docs entries. With the certified kernels the locked steady state is well posed (the 1061 x 1061 block has condition number 3e5, all eigenvalues negative; the Tikhonov weight changed the result by at most 0.1% from 0.01 to 10), and it is the steady point of the semi-analytical system, `y = -A \ c`, so no separate solver is needed.*
- [x] ~~`filenameGenerator.m`: no callers in the tracked code.~~ *Kept on purpose:* it gives each run its own numbered folder (`000001/results.mat`, `000002/…`), which suits batch jobs of many forward models (strain rates, stresses, surface observations) whose outputs should be stored systematically. `experiments/EXPT1_run.m` uses it this way.
- [ ] Optional tidy-up of `filenameGenerator.m`: the long `if/elseif` chain could be one line, `sprintf('%06d', num)`.
- [x] ~~`pcolor_logvals.m`: see §1. *Deleted.*~~

## 4. File and function naming

**Rules (adopted 2026-11-21).**
1. Functions: lowerCamelCase, starting with a verb: `compute…`, `assemble…`, `solve…`, `run…`, `integrate…`, `extract…`, `plot…`, `load…`/`save…`, `create…`, `certify…`, `is…`. No snake_case, no `get…`.
2. Acronyms are one capitalised word (`Bem`, `Shz`, `Ivp`, `Ode`), except long-established ones (`LD…`, `KK`/`KL`/`LK`/`LL`).
3. Half-space kernels live in `+geometry/` with no suffix; the boundary-corrected versions live in `functions/` with the suffix `Bem` (e.g. `geometry.computeFaultTractionKernels` vs `computeFaultTractionKernelsBem`).
4. Domain words: `Fault` (megathrust patches), `Shz` (mantle cells, as the `shz` variable), `Boundary` (model edge).
5. Classes keep their lowercase names (`geometry.receiver`, `shearZoneReceiver`, `earthModel`, `fpatch`, `shzpatch`, `LDhs`): every script uses them, and on a case-insensitive file system `receiver.m` and `Receiver.m` cannot coexist, so there is no way to keep the old name working. A documented exception.
6. Third-party files keep their names: `LDdispHS`, `LDstressHS`, `cline`, `bluewhitered`, `drchrnd`. `filenameGenerator` is kept as is (used by experiment scripts).
7. Scripts: `exampleN_whatItShows`, `plot…`.

**Renamed (2026-11-21).** Every caller, test, example and docs page was updated; scripts outside the repository (e.g. the git-ignored `experiments/`) were not.

| Old | New |
| --- | --- |
| `computePeriodicInitialConditionBalanced` | `computePeriodicInitialCondition` |
| `eigBalanced` | `eigOperator` |
| `analyticSolveIVP` | `solveAnalyticIvp` |
| `egh_modal` | `computeModalIntegrals` |
| `ode_viscoelastic` | `odeViscoelastic` |
| `viscoeq_sequence_simulation` | `integrateImposedCycles` |
| `extract_v_edot_from_statevector` | `extractRatesFromOdeState` |
| `extractComponentsSolutionVector` | `extractRatesFromLinearState` |
| `getStrainratesLongterm` | `computeLongTermStrainRates` |
| `computeCornerflowVelocityfield` | `computeCornerFlowVelocity` |
| `cornerFlowVelocity` | `computeRegionCornerFlow` |
| `plotpatch2d` | `plotPatches2d` |
| `plotshz2d` | `plotShz2d` |
| `certifyMantleKernel` | `certifyShzKernel` |
| `geometry.flt2flt` | `geometry.splitSegment` |
| `create_megathrust` | `createMegathrust` |

**Deleted (§3):** `computePeriodicInitialCondition` (the legacy expm/SVD version), `reconditionODEoperator`, `computeKLmodesMesh`, `computeSiteVelocitiesBem`.

Error identifiers follow their functions (`certifyShzKernel:notNegative`, `certifyShzKernel:zeroEnergyModes`, `computeRegionCornerFlow:region`).

- [x] ~~Credit `bluewhitered`, `cline` and `drchrnd` under "Third-party components" in `docs/contributors.md`.~~ *Done: bluewhitered (Nathan Childress, MATLAB File Exchange 4058, 2003) and cline (Sebastian Hölz, File Exchange 14677, 2007) are credited with links; drchrnd is described as the standard gamma-variate construction. The naming rules are in the docs (Contributors → Code conventions) and the old → new table is on the "Renamed functions" page.*

## 5. `meshing/` folder

- [x] ~~Replaced the vendored `mesh2d` library and the two outdated scripts (`CREATE_shearzone_mesh.m`, `CREATE_geometryFromEdgesMesh.m`; neither produced the mesh in `inputs/`, which covers a larger domain) with one function, `meshing/createShearZoneMesh.m`. It meshes the example outline with the PDE Toolbox and writes the `.dat` files. Documented in `docs/inputs.md`; test: `GeometryTest/createdMeshIsReadableAndFillsTheOutline`.~~
- [x] ~~The exact settings of the example mesh are lost (`'Hmax', 20` gives 495 triangles instead of 521).~~ *Decided: keep the mesh in `inputs/` as it is; `createShearZoneMesh` is only for making new meshes.*

## 6. Cached kernels and data

- [ ] `examples/kernels/evl_orig.h5` (~8.7 MB) and `examples/kernels/devl.h5` (~6.2 MB) are loaded by `examples/exampleModel.m` (for examples 1, 2 and 4) and `tests/test_negate_exact_solver.m`, and nothing checks that they still match the current `examples/inputs/` mesh or kernel code beyond their sizes. → Add a regeneration script, plus a check on mesh size or a hash. Or regenerate them when they're missing and stop tracking them (Git LFS or a release asset).
- [x] ~~Store kernels in a format usable outside MATLAB.~~ *Done (2026-11-20): `saveKernels`/`loadKernels` write and read HDF5 (format 'I2EC kernels', version 1). Each kernel block is a compressed dataset with description and units attributes; the certificate, the build options (`evl.options`, now recorded by `computeAllStressKernelsBem`) and optionally the fault and mantle geometry are stored too. Row-major readers (h5py/NumPy) see the MATLAB shape and index order. `kernels/evl_orig.mat` and `devl.mat` were converted to `.h5` (identical on reading back) and removed; example 3 caches `.h5` and converts an existing `.mat` cache once. Layout: `docs/inputs.md`. Tests: `KernelsTest/kernelFilesRoundTrip`, `kernelFilesShowTheMatlabShapeToRowMajorReaders`, `loadKernelsRejectsOtherFiles`.*
- [x] ~~`MAIN_analyticsolve_singleeqcycle.m` keeps the `computeAll*KernelsBem` calls as comments and loads the `.mat` files instead. → Replace this with an explicit `if exist(...) load else compute` switch.~~ *Done: `examples/exampleModel.m` loads the cached kernels when they exist and match the mesh sizes, and computes them otherwise.*

## 7. Path handling and project setup

- [ ] Every script runs its own `addpath functions/` and depends on the current folder (`load('kernels/...')`). → Add one `setup.m` or `startup.m` that adds `functions/` and the repo root (for `+geometry`) with `fullfile(fileparts(mfilename('fullpath')), ...)`, and call it from the scripts and the tests.
- [x] ~~`PLOT_longtermVelocity.m` adds `inputs/` to the path. Data folders shouldn't be on the MATLAB path; load files by full path instead.~~ *Done: `examples/plotLongtermVelocity.m` loads its inputs by full path through `exampleModel`.*

## 8. Tests and CI

- [x] ~~Test suite in `tests/` (`matlab.unittest`, runner `tests/runTests.m`, guide `tests/README.md`): 84 fast tests (Green's functions, geometry, BEM kernels, corner flow, solvers, ODE path, utilities, legacy code), plus a slow tier that rebuilds the example kernels and runs the MAIN scripts.~~
- [ ] Run `runTests` in CI with `matlab-actions/setup-matlab` + `run-command`. The runner exits non-zero only on unexpected failures, so `KnownIssue` tests don't block CI. Needs a MATLAB licence in CI (free for public repos via MathWorks' GitHub actions); Parallel Computing Toolbox is optional (`parfor` runs serially without it).
- [x] ~~Remove each `KnownIssue` tag as its bug in §0 is fixed.~~ *Done: no `KnownIssue` tests remain.*
- [ ] The docs workflow (`.github/workflows/docs.yml`) only builds Sphinx (`docs/build.sh`). `docs/_ext/check_help.py` exists and exits 1 on a failed check, but only runs by hand. → Run it in CI so every documented function exists and has a help block.

## 9. Documentation clean-up (after the deletions above)

- [x] ~~Remove the docs entries for deleted functions in `docs/api/kernels.rst`, `solvers.rst`, `utilities.rst` and `longterm.rst`.~~ *Every `mat:autofunction` entry now points at an existing function (checked 2026-11-20).*
- [x] ~~Bring help blocks into one format (H1 line, INPUTS/OUTPUTS, author). Several files have no H1 line (for example `computeSiteVelocitiesBem`, `assembleStressKernel`), or have a stale name in it (`computeKLmodesMesh`).~~ *Done (checked 2026-11-20): every file in `functions/` and `+geometry/` has an H1 line naming itself, and every file in `functions/` has INPUTS/OUTPUTS.*
- [x] ~~Bring the README's "Code Dependencies" list up to date (bluewhitered, cline, drchrnd sources).~~ *Superseded: the README was rewritten around the documentation website (#13) and has no dependency list. Credits live in `docs/contributors.md`; the missing ones are tracked in §4.*

## 11. Analytic stress of a strained triangle (2026-11-20, merged into main in #14)

**Problem.** `computeStressPlaneStrainTriangleShearZoneFiniteDifference` gets the strain by differencing the closed-form displacement (step 1e-6 l, one-sided in x3): four displacement evaluations per point, truncation errors of ~1e-6, and a patch of special cases (`atan(N/D)` jumps by pi on lines through source and image edges, d5279e7). The displacement is a closed form, so its gradient is one too.

**What was done.**
- `derivation/generateTriangleStress.m` differentiates the edge integrals `I22`, `I23`, `I32`, `I33` (every branch: general, a2 = 0, a3 = 0; read from the displacement file by line number) with the Symbolic Math Toolbox and writes `functions/computeStressPlaneStrainTriangleShearZone.m` (hand-written part: `derivation/triangleStressTemplate.m`; generated local functions below it, ~1 min to regenerate). d/dx atan(N/D) is written (N'D - ND')/(N^2 + D^2).
- Same signature as the FD routine; e22, e23, e33 may also be [1 x K] rows (K strains share the edge integrals).
- **Removable singularities (found in testing).** Terms like atan(N/D)/D cancel only between the two ends of an edge. On the line through an *image* edge, and through a *horizontal or vertical* source edge (the a2 = 0 / a3 = 0 branches), the closed form is NaN, and nearby it loses accuracy like 1e-13*(r/d)^2 (d: distance to the line, r: distance from the edge; measured 1e6 relative error at d = 1e-10 L, 1e-5 at 1e-4 L, and 3e-6 at d = 0.02 L when r = 45 L). General source-edge lines are fine. Fix: for points with d < 1e-2 max(L, r), that edge's contribution is integrated numerically (differentiated Green's function along the edge: 64-point Gauss-Legendre farther than L/4 from the edge, tanh-sinh nearer). Non-finite results (on the triangle's own edges/vertices) also get quadrature.
- The Green's-function derivative uses atan2(x2-y2, x3+y3), so it stays finite at x3 = y3 = 0 (edge on the surface, point on the surface).
- `geometry.computeShzStressKernels(src, shz)` uses it, one call per source for all three strains. (It briefly took a `'finitedifference'` option; removed with the FD routine, see Left.)
- Docs: `docs/theory.md` ("Mantle stress from the differentiated displacement"), `docs/api/kernels.rst`.

**Results.**
- [x] ~~Independent reference (tanh-sinh quadrature, step 2e-4, of the differentiated Green's function; scratch script, not in the repo)~~: 16 random triangles incl. axis-aligned edges, ~280 points each packed onto every source and image edge line (d from 1e-8 to 1e-2 L) and random points: median 1e-12, worst 2e-9 relative. FD against the same reference: median 1e-6, max 1.4e-2 (a shallow triangle).
- [x] ~~Traction-free surface to ~1e-14 (also with an edge on the surface).~~ Equilibrium (Richardson FD of the stress, |div s| l/|s|): analytic median 1e-9, FD 1e-4. Betti reciprocity with exact cell averages (16x16 Gauss): 1e-14 to 6e-12 analytic, 1e-6 FD, at 2.5-30 km separation.
- [x] ~~Smooth across the edge lines without the FD routine's shifted-point special case~~; continuous to ~1e-8 at 1e-6 m offsets.
- [x] ~~Example mesh, half-space `LL` (3-point)~~: entries agree with FD to 1.5e-8 median, 8e-5 max, of each source's largest stress; the worst disagreements were checked against the quadrature reference and the analytic value is right in every case (e.g. FD -2.49 vs analytic and reference -0.0147). 0 growing modes with both (max Re/max abs -2.2e-3).
- [x] ~~Full example kernel (`computeAllStressKernelsBem`, certified 3-point) vs cached FD `evl_orig.h5`~~: KK, KL identical; LK 7e-6, LL 4e-5 of max, within `cachedStressKernelsAreReproducible`'s then 1e-4 (the cache was since regenerated and the test tightened to 1e-8, see Left). Certifies with 0 clipped modes, max eigenvalue -8.3e-4 of the spectrum.
- [x] ~~Aligned shallow structured mesh (diagonal, 15 km, z -150 to -45 km, 280 cells; the d5279e7 failure case)~~: analytic = FD after that fix: centroid 43 growing modes (max Re/max abs 0.14), 3-point 0.
- [x] ~~Speed~~: per source triangle with 1563 receivers, 928 ms FD vs 85 ms analytic (11x). Half-space shz->shz kernel: 23.7 s → 2.5 s. `computeShzStressKernelsBem` was then dominated by receiver-independent work (source→boundary kernels and one BEM solve, plus a full-matrix product, per source), repeated for each of the 3 quadrature points and again for LK; fixed below.
- [x] ~~Tests~~: 6 new tests in `GreensFunctionsTest` (vs a 4th-order derivative of the displacement, vs FD, traction-free surface, Betti, smooth on edge lines, several strains at once). Fast tier: 141 passed, 0 unexpected failures (3 known issues). The vs-FD test went with the FD routine; the derivative test now also checks a point inside the triangle. Fast tier after the merge: 140 passed.

**Left.**
- [x] ~~Regenerate the kernel caches with the analytic stress.~~ *Done after merging `maintenance-cleanup` (HDF5). `kernels/evl_orig.h5`: KK, KL unchanged, LK 7e-6, LL 4e-5 of max; certified, 0 clipped. `kernels/devl.h5` also rebuilt (its mantle part uses the stress for the traction BC): KO 6e-12, LO 6e-10. Mesh-convergence caches (git-ignored): h40 / h20 / h10 in 12 / 13 / 51 s (the fine mesh took about an hour before), LL change 5e-6 / 4e-5 / 6.5e-4 (FD error grows as cells shrink), all certified, 0 clipped. `KernelsTest/cachedStressKernelsAreReproducible` now holds LK, LL to 1e-8 like KK, KL; passes. Fast tier 143 passed, 0 unexpected failures.*
- [x] ~~Merge the branch into main (PR #14).~~ *Merged 2026-11-20, with its docs revision (#15).*
- [x] ~~Decide on the FD routine.~~ *Deleted (2026-11-20): `computeStressPlaneStrainTriangleShearZoneFiniteDifference` and the `'finitedifference'` option of `geometry.computeShzStressKernels`. Its tests went too (traction-free surface, Betti, analytic-vs-FD), each covered by an analytic counterpart; the FD edge-line regression case now runs on the analytic routine, and the displacement-derivative test gained a point inside the triangle (the only interior check the analytic-vs-FD test gave).*
- [x] ~~Follow-up (speed): compute the receiver-independent BEM part once.~~ *Done: new `computeShzBoundaryCorrection(src, boundary)` returns the boundary slip that cancels each unit source strain (`.slip`, [2Nb x Ns x 3]; one LU, all 3Ns right-hand sides). `computeShzStressKernelsBem(src, shz, boundary, bc)` takes it as an optional 4th argument (computed inside if omitted; size-checked) and is now matrix products only. `computeAllStressKernelsBem` computes it once and shares it between the 3 quadrature points and LK. Example mesh: shz→shz 203 s → 10 s (7 s with `bc` passed, `bc` itself 6 s), shz→fault 193 s → 1.4 s; results equal to the old loop to 1e-17. Full `computeAllStressKernelsBem` (certified 3-point): 27 s, equal to the earlier analytic build to 2e-17. Tests: `BoundaryCorrectionTest/precomputedCorrectionGivesSameKernels`, `correctionForOtherSourcesIsRejected`. Fast tier 140 passed, 0 unexpected failures. The fault-source routines (`computeFaultTractionKernelsBem`, KL block) still solve per source; not measured as a bottleneck.*
- [ ] Optional: the generator reads the displacement file by line number; regenerate after any edit to `computeDisplacementPlaneStrainTriangleShearZone`.

---

### Suggested order

0. §0: fix the bugs in code that is in use (`getStrainratesLongterm` done).
0a. §0a: done; certified 3-point kernels are the default (`kernelmodify` kept as the earlier fix).
1. §1 and §2: delete or fix the broken files and duplicates (low risk, since nothing calls them).
2. §3 and §5: decide which unused files to keep.
3. §7: one setup script.
4. §4: renames, done after the code base is smaller.
5. §6, §8, §9: data, tests, docs.
