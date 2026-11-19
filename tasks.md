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

---

## 0. Bugs confirmed by the test suite (2026-11-19)

Each item names the failing test (`Class/test`).

### Affects code that is in use

- [x] ~~**`MAIN_i2ec_singleeqcycle.m` and `computeInterseismicSteadystate.m` cannot run.** *Fixed: `n` defaults to 1, and both callers now pass `n = 1`.* Both call `getStrainratesLongterm(shz, dip, oc, co)` with 4 arguments, but the function now requires the exponent `n`, so they fail with "Not enough input arguments". → Give `n` a default of 1 in `getStrainratesLongterm`, or pass it at both call sites. Tests: `LongTermFlowTest/strainRatesWithoutExponentDefaultToNewtonian`, `SolversTest/interseismicSteadyStateRuns`, `ExampleScriptsTest/numericalSingleCycleScriptRuns`.~~
- [x] ~~**`runImposedEarthquakeCycles(..., nreps = 1, ...)` fails.** The first call discards `t` (`[\~,Y] = ...`) and the loop never runs, so output `t` is never assigned. → Keep `t` from the first call. *Fixed: the first call now keeps `t`, so `nreps = 1` returns one cycle from the long-term state. `nreps > 1` is unchanged. Test: `TimeIntegrationTest/singleCycleRunReturnsTime` (no longer `KnownIssue`; also checks `t(end) = Trecur` and the output sizes). Since then `viscoeq_sequence_simulation` takes `nreps` and runs every segment of every cycle in one loop, so `runImposedEarthquakeCycles` no longer has a cycle loop of its own. Results are bitwise identical to the old double loop (checked for 1–7 events, 1–3 cycles, n = 1 and 3, with locked patches). Tests: `TimeIntegrationTest/cyclesChainFromTheEndState`, `everyEventIsAppliedOnce`, `badEventTimesAreRejected`.*~~
- [ ] **The megathrust mesh has one patch fewer than asked for, and the deepest patch is twice as wide.** `flt2flt` merges the last two patches when the remaining length is exactly `2*W0`, because the floating-point `Wt > Wc/2` comparison is true for equal values. So `megathrust2d.seg` (W0 = W/50) gives 49 patches, and `create_megathrust(..., Nfault)` gives `Nfault-1`. → Use a tolerance in the comparison, or compute the patch count first. Rebuild `kernels/*.mat` afterwards (§6). Test: `GeometryTest/uniformSegmentGivesUniformPatches`.
- [x] ~~**`geometry.computeFaultDisplacementKernels` nudges a point that sits on a patch centre along the wrong normal.** It uses `rcv.nv(r==0,:)`, where `r==0` is a mask over observation points, not patches. The point gets another patch's normal, or the call errors ("true value outside of the array bounds") when the point's row number is greater than `rcv.N`. → Use `rcv.nv(k,:)`. Also, `delta = 1e-9*mean(ox(:,1))` is 0 or negative when the mean x is ≤ 0. *Fixed: a point on the centre of patch k now moves by `1e-9*rcv.W(k)` along `rcv.nv(k,:)`, so it always gets the value on the +nv side. No current result changes. The only callers that put points on patch centres are the BEM self-problems (`boundary` on `boundary.xc`), where the old indexing happened to pick the right normal and the step was positive. Every BEM kernel of the example (`evl`, `devl`, the hinge kernels) matches the old values to 1e-8 or better, so `kernels/*.mat` don't need regenerating. Tests: `GeometryTest/faultDisplacementKernelOnPatchCentre` (now also checks the value) and `faultDisplacementKernelOnPatchCentreNegativeX` (new). Both fail on the old code.*~~
- [ ] **The BEM stress kernels depend on the boundary velocity.** `computeFaultTractionKernelsBem` and `computeShzStressKernelsBem` put `boundary.Vx`/`Vz` into the Dirichlet data of every unit-source column. The displacement-kernel routines use `0.*Vx`. With `bcx`/`bcz ≠ 0` in the boundary `.seg` file (for example `boundarytight2d.seg`), every kernel column picks up a constant term, so the kernels are no longer linear in the source. This has no effect on the current example, because `boundary2d.seg` has zero boundary velocity. → Use `0.*Vx` there too, and treat a prescribed boundary motion as a separate loading term. Test: `KernelsTest/stressKernelsIgnoreBoundaryVelocity`.
- [x] ~~**The n = 3 corner flow runs `fsolve` on every call, and prints its output each time.** It also never checked convergence: for dips below about 8° the wedge solve failed silently and returned velocities off by up to 0.6 of the plate rate. *Fixed: the four coefficients now come from one scalar equation for the phase (`fzero`, every root bracketed), with the other three solved linearly; the result is cached per dip and an error is raised if no root meets the boundary conditions. No Optimization Toolbox needed. Tests: `LongTermFlowTest/powerLawCornerFlowBoundaryConditions` and `cornerFlowSolvesStokes` over dips 3–80°, `stokesCheckRejectsWrongRheology`, `powerLawCornerFlowIsQuiet`, `powerLawCoefficientsAreSolvedOncePerDip`.*~~
- [x] ~~**`computeShzStressKernelsBem` added the box correction to the wrong stress components.** For mantle receivers, the half-space stress is stacked as (σxx, σxz, σzz), but the BEM box correction from `geometry.computeFaultStressKernels` was stacked as (σxx, σzz, σxz). So the σxz and σzz corrections were swapped in `LL`. The full-kernel reciprocity test missed it because the correction is only \~1% of `LL`; the correction alone was 65% asymmetric under Betti's theorem. *Fixed in `2eede92` by reordering the rows, and `kernels/evl_orig.mat` was regenerated. `LL` changes by 0.5%, and `KK`, `KL` and `LK` are unchanged. The growing modes and the fastest decay rates don't move, but the slow modes do: in `MAIN_analyticsolve_singleeqcycle.m`, strain rates 50 years after the earthquake change by up to \~20%, and late-interseismic surface velocities by up to 0.04 Vpl. Results computed before the fix carry this error. Test: `BoundaryCorrectionTest/boxCorrectionIsReciprocal` (correction alone, 0.15% asymmetric after the fix).*~~

### Affects unused or legacy code (delete rather than fix, see §1-§3)

- [x] ~~**`computeStressPlaneStrainTriangleShearZone.m` (analytic, 2,238 lines) gives wrong stresses.** It disagrees with the finite-difference version everywhere, including the traction-free surface. Line 84 overwrites the summed `u2d2` with a single term, which looks like leftover debugging. The eigenstrain is then removed using `e22`/`e23`/`e33` after they have been overwritten with the total strain. And `Omega(x2,x3)*e22` is a matrix product, so vector input errors. → Delete it. The finite-difference version is correct: it passes the surface-traction and reciprocity tests. Tests: `GreensFunctionsTest/triangleStressAnalytic*`. *Deleted, with its two tests and its docs entry. Nothing ever called it.*~~
- [x] ~~`computeStressFault_bem.m` errors with "Undefined function `computeFullTractionKernels`". *Deleted, with its docs entry and its test.* Test: `LegacyCodeTest/stressFaultBemRuns`.~~
- [x] ~~`getStrainratesLongterm_old.m` declares `getStrainratesLongterm`. *Deleted.* Test: `LegacyCodeTest/functionNamesMatchFileNames`.~~
- [x] ~~Four file names exist in both `functions/` and `+geometry/`. *Fixed: the `functions/` copies are deleted.* Test: `LegacyCodeTest/noDuplicateFileNames`.~~
- [ ] `reconditionODEoperator.m`: the help text says growing eigenvalues are replaced by `-eps`, but the code uses `-0.1` s⁻¹, a 10-second relaxation. `computePeriodicInitialCondition.m` prints the condition number on every call. Tests: `SolversTest/reconditionReplacesGrowingModesByMinusEps`, `SolversTest/legacyPeriodicSolverIsQuiet`.
- [ ] `pcolor_logvals.m` draws nothing. Test: `UtilitiesTest/pcolorLogValsDrawsSomething`.
- [ ] `create_megathrust.m` writes `megathrust2d.seg` into the current folder. Test: `UtilitiesTest/createMegathrustLeavesWorkingFolderClean`.

### Discretisation (known, handled by the solvers)

- [ ] The collocated mantle-mantle block `LL` has eigenvalues with positive real part: 43 on the example mesh (max 2.4e3), and some even on a 24-triangle mesh. `negateGrowingModes` and `kernelmodify = 1` handle them. A better quadrature for the self-term and nearest neighbours might remove them at the source. Test: `KernelsTest/unmodifiedMantleKernelHasNoGrowingModes`.

### Confirmed correct (safe to rely on during the clean-up)

- The kernel sign conventions agree with each other: KK, KL/LK and LL satisfy Betti reciprocity (ratio ≈ +1) on both the example and the synthetic mesh. This full-kernel check is dominated by the half-space part. The BEM box correction is now checked separately by `BoundaryCorrectionTest`, which caught the swapped rows above. `computeFaultStressKernels`' `-Sxz` convention matches `computeFaultTractionKernels`.
- **Old corner flow = new corner flow** (n = 1) to 1e-14. So `compute_cornerflow_velocityfield` can be replaced by `computeCornerflowVelocityfield(…, 1)` in `computeLongTermVelocity` without changing results (item in §1). Both satisfy their boundary conditions and are divergence-free.
- **Duplicate copies give identical results**: `functions/LDdispHS`, `EdgeStress`, `computedisplacementGFs`, `computetractionkernels`, and the old edge-dislocation traction kernels, which match the line-dislocation kernels. The deletions in §2 change no numbers.
- `computeSiteVelocitiesBem` gives the same result as `devl` from `computeAllDisplacementKernelsBem`, to 1e-12.
- **The ODE path (`ode_viscoelastic` + `ode45`) matches `analyticSolveIVP`** for a locked fault with a Newtonian mantle, through a coseismic step. The solvers match `expm` to round-off.

## 1. Broken or likely incorrect code (fix or delete first)

- [x] ~~**`functions/computeStressFault_bem.m` is broken.** *Deleted, with its docs entry and `LegacyCodeTest/stressFaultBemRuns`.* It called `computeFullTractionKernels`, which doesn't exist anywhere in the repo. Nothing called this function, but the docs listed it (`docs/api/kernels.rst`).~~
- [x] ~~**`functions/getStrainratesLongterm_old.m` has the wrong function name inside it.** *Deleted.* The file declares `function ... = getStrainratesLongterm(...)`, so the declared name doesn't match the filename. Nothing calls it. → Delete it.~~
- [ ] **`functions/pcolor_logvals.m` is unfinished.** It builds `pos_logvals` and `neg_logvals` and then returns without plotting or outputting anything. Nothing calls it, but the docs list it (`docs/api/utilities.rst`). → Finish it or delete it.
- [ ] **Two corner-flow implementations.** `compute_cornerflow_velocityfield.m` (old, JPL 2023, Newtonian only) is still used by `computeLongTermVelocity.m` (and so by `PLOT_longtermVelocity.m`). `getStrainratesLongterm.m` uses the newer `computeCornerflowVelocityfield.m` (supports n = 1 and n = 3). The test suite shows the two agree to 1e-14 for n = 1 (`LongTermFlowTest/oldAndNewNewtonianCornerFlowAgree`). → Switch `computeLongTermVelocity` to `computeCornerflowVelocityfield(…, 1)` and delete the old one.
- [ ] **`computeLongTermVelocity.m` and `getStrainratesLongterm.m` repeat the same translate/rescale block** (`ocean_translate_params`, `continental_translate_params`, the `- 1` shift). → Pull it into one helper so the two can't drift apart.
- [ ] **`computeInterseismicSteadystate.m` hard-codes geometry inputs:** `[0,20e3]` and `[-140e3,35e3]` are passed to `getStrainratesLongterm`. → Make them arguments, or document them. Moot if the function is deleted (§3): since the notebook was removed, nothing in the tracked code calls it.

## 2. Duplicate code: `+geometry/` vs `functions/`

Several files exist in both places. Package code calls the `geometry.` versions, so the `functions/` copies are leftovers from older versions.

| File | `+geometry/` copy | `functions/` copy | Action |
|---|---|---|---|
| `LDdispHS.m` | used (`geometry.LDdispHS`) | ~~identical, unused~~ **deleted** | done |
| `EdgeStress.m` | only used by the unused `geometry.computetractionkernels` | ~~identical, unused~~ **deleted** | delete the `+geometry/` copy together with `geometry.computetractionkernels` |
| `computedisplacementGFs.m` | **unused** | ~~unused~~ **deleted** | delete the `+geometry/` copy (replaced by `computeFaultDisplacementKernels`) |
| `computetractionkernels.m` | **unused** | ~~unused~~ **deleted** | delete the `+geometry/` copy (replaced by `computeFaultTractionKernels`) |

- [x] ~~Delete the `functions/` copies of the four files above.~~
- [ ] Delete the unused `+geometry/` versions: `computedisplacementGFs`, `computetractionkernels` and `EdgeStress` (only `computetractionkernels` uses it). `LegacyCodeTest` shows they match the current kernels, so no numbers change. Delete their equivalence tests in `LegacyCodeTest` along with them. Their docs entries go with the files, because `docs/api/geometry.rst` documents the package automatically.
- [x] ~~Remove the docs mentions of the `functions/` copies (the note in `docs/api/kernels.rst`; they had no API entries of their own).~~
- [ ] Decide on one rule for what goes where. For example: low-level Green's functions and geometry classes go in `+geometry/`, and the `*Bem` assembly, solvers and plotting go in `functions/`. Write the rule down in the README or docs.

## 3. Unused functions in `functions/`

No caller in the repo. Each one needs a keep-or-delete decision.

- [ ] `computeKLmodesMesh.m`: KL/eigenmode basis, looks experimental. Also, the help text calls it `compute_eigmodes_mesh` (out-of-date name).
- [ ] `computePeriodicInitialCondition.m` (with `reconditionODEoperator.m`, which only it uses): replaced by `computePeriodicInitialConditionBalanced.m` and `negateGrowingModes.m`. → Delete both, or keep one as a documented fallback.
- [ ] `computeSiteVelocitiesBem.m`: nothing calls it, and it carries a block of commented-out `devl`-based code.
- [x] ~~`computeStressPlaneStrainTriangleShearZone.m` (**2,238 lines**, Barbot 2018 analytic stress): the code actually uses `computeStressPlaneStrainTriangleShearZoneFiniteDifference.m` instead. The tests show it gives **wrong stresses** (§0), so it cannot serve as a reference. → Delete it. This is the largest single file in the repo. *Deleted.*~~
- [ ] `computeInterseismicSteadystate.m`: its only caller was the notebook, now deleted (the remaining callers are in the git-ignored `experiments/`). Interseismic steady state with a pinned fault and a Laplacian-smoothed mantle; the semi-analytical driver does not use it. → Delete, or keep and add it to a driver script.
- [ ] `filenameGenerator.m`: no callers. Its long `if/elseif` chain could be one line, `sprintf('%06d', num)`.
- [ ] `pcolor_logvals.m`: see §1.

## 4. File and function naming

The repo mixes several naming styles:

- **camelCase**: `computeLongTermVelocity`, `analyticSolveIVP`, `negateGrowingModes`
- **snake_case**: `compute_cornerflow_velocityfield`, `extract_v_edot_from_statevector`, `create_megathrust`, `viscoeq_sequence_simulation`, `ode_viscoelastic`, `egh_modal`, `pcolor_logvals`
- **all lowercase**: `computedisplacementGFs`, `computetractionkernels`, `bluewhitered`, `drchrnd`, `cline`, `plotpatch2d`, `plotshz2d`, `fpatch`, `shzpatch`, `flt2flt`, `receiver`
- **UPPERCASE script prefixes**: `MAIN_`, `PLOT_`, `CREATE_`, `TEST_`, plus a lowercase `tests/test_*.m`

Tasks:

- [ ] Pick a convention and write it in the README or docs. Suggestion: lowerCamelCase for functions, UpperCamelCase for classes (`Receiver`, `ShearZoneReceiver`, `EarthModel`), and `run*` / `plot*` or a `scripts/` or `examples/` folder for entry-point scripts.
- [ ] Rename the snake_case project functions after the dead code is gone, so you don't rename files that are about to be deleted. Do it in one commit per group, and update every caller and `docs/api/*.rst`. Candidates:
  - `extract_v_edot_from_statevector` → `extractVEdotFromStateVector`
  - `viscoeq_sequence_simulation` → `simulateViscoEqSequence`
  - `ode_viscoelastic` → `odeViscoelastic`
  - `create_megathrust` → `createMegathrust`
  - `egh_modal` → `eghModal` (or a more descriptive name)
- [ ] Decide whether to rename the `+geometry` classes (`receiver`, `shearZoneReceiver`, `earthModel`, `fpatch`, `shzpatch`). Renaming them changes the public API for anyone with existing scripts.
- [ ] Leave third-party files under their original names: `bluewhitered`, `cline`, `drchrnd`, `LDdispHS`, `LDstressHS`, `EdgeStress`. Credit their source in a header or a `THIRD_PARTY.md`.
- [ ] Bring the example and entry-point scripts in line with each other (`MAIN_analyticsolve_singleeqcycle.m`, `MAIN_i2ec_singleeqcycle.m`, `PLOT_longtermVelocity.m`). Consider moving them into `examples/`.

## 5. `meshing/` folder

- [x] ~~Replaced the vendored `mesh2d` library and the two outdated scripts (`CREATE_shearzone_mesh.m`, `CREATE_geometryFromEdgesMesh.m`; neither produced the mesh in `inputs/`, which covers a larger domain) with one function, `meshing/createShearZoneMesh.m`. It meshes the example outline with the PDE Toolbox and writes the `.dat` files. Documented in `docs/inputs.md`; test: `GeometryTest/createdMeshIsReadableAndFillsTheOutline`.~~
- [ ] The exact settings of the example mesh are lost: `'Hmax', 20` gives 495 triangles instead of 521. → Either keep `inputs/` as is, or regenerate it with the new function and rebuild `kernels/` (changes the example's numbers).

## 6. Cached kernels and data

- [ ] `kernels/evl_orig.mat` (~9 MB) and `kernels/devl.mat` (~6.7 MB) were committed in Jan 2024. They're loaded by `MAIN_analyticsolve_singleeqcycle.m` and `tests/test_negate_exact_solver.m`, and nothing checks that they still match the current `inputs/` mesh or kernel code. → Add a regeneration script, plus a check on mesh size or a hash. Or regenerate them when they're missing and stop tracking them (Git LFS or a release asset).
- [ ] `MAIN_analyticsolve_singleeqcycle.m` keeps the `computeAll*KernelsBem` calls as comments and loads the `.mat` files instead. → Replace this with an explicit `if exist(...) load else compute` switch.

## 7. Path handling and project setup

- [ ] Every script runs its own `addpath functions/` and depends on the current folder (`load('kernels/...')`). → Add one `setup.m` or `startup.m` that adds `functions/` and the repo root (for `+geometry`) with `fullfile(fileparts(mfilename('fullpath')), ...)`, and call it from the scripts and the tests.
- [ ] `PLOT_longtermVelocity.m` adds `inputs/` to the path. Data folders shouldn't be on the MATLAB path; load files by full path instead.

## 8. Tests and CI

- [x] ~~Test suite in `tests/` (`matlab.unittest`, runner `tests/runTests.m`, guide `tests/README.md`): 84 fast tests (Green's functions, geometry, BEM kernels, corner flow, solvers, ODE path, utilities, legacy code), plus a slow tier that rebuilds the example kernels and runs the MAIN scripts.~~
- [ ] Run `runTests` in CI with `matlab-actions/setup-matlab` + `run-command`. The runner exits non-zero only on unexpected failures, so `KnownIssue` tests don't block CI. Needs a MATLAB licence in CI (free for public repos via MathWorks' GitHub actions); Parallel Computing Toolbox is optional (`parfor` runs serially without it).
- [ ] Remove each `KnownIssue` tag as its bug in §0 is fixed.
- [ ] The docs workflow (`.github/workflows/docs.yml`) only builds Sphinx. → Add `docs/_ext/check_help.py`, or a similar check, to CI so every documented function exists and has a help block.

## 9. Documentation clean-up (after the deletions above)

- [ ] Remove the docs entries for deleted functions in `docs/api/kernels.rst`, `solvers.rst`, `utilities.rst` and `longterm.rst`.
- [ ] Bring help blocks into one format (H1 line, INPUTS/OUTPUTS, author). Several files have no H1 line (for example `computeSiteVelocitiesBem`, `assembleStressKernel`), or have a stale name in it (`computeKLmodesMesh`).
- [ ] Bring the README's "Code Dependencies" list up to date (bluewhitered, cline, drchrnd sources).

---

### Suggested order

0. §0: fix the bugs in code that is in use (`getStrainratesLongterm` done).
1. §1 and §2: delete or fix the broken files and duplicates (low risk, since nothing calls them).
2. §3 and §5: decide which unused files to keep.
3. §7: one setup script.
4. §4: renames, done after the code base is smaller.
5. §6, §8, §9: data, tests, docs.
