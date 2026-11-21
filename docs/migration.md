# Renamed functions

In November 2026 sixteen functions were renamed to one naming scheme (see
{ref}`code-conventions`) and four unused ones were deleted. The old names no
longer exist, so scripts written for earlier versions need updating.

## Renamed

Each function keeps its arguments and outputs; only the name changed.

| Old name | New name |
| --- | --- |
| `computePeriodicInitialConditionBalanced` | {mat:func}`computePeriodicInitialCondition` |
| `eigBalanced` | {mat:func}`eigOperator` |
| `analyticSolveIVP` | {mat:func}`solveAnalyticIvp` |
| `egh_modal` | {mat:func}`computeModalIntegrals` |
| `ode_viscoelastic` | {mat:func}`odeViscoelastic` |
| `viscoeq_sequence_simulation` | {mat:func}`integrateImposedCycles` |
| `extract_v_edot_from_statevector` | {mat:func}`extractRatesFromOdeState` |
| `extractComponentsSolutionVector` | {mat:func}`extractRatesFromLinearState` |
| `getStrainratesLongterm` | {mat:func}`computeLongTermStrainRates` |
| `computeCornerflowVelocityfield` | {mat:func}`computeCornerFlowVelocity` |
| `cornerFlowVelocity` | {mat:func}`computeRegionCornerFlow` |
| `plotpatch2d` | {mat:func}`plotPatches2d` |
| `plotshz2d` | {mat:func}`plotShz2d` |
| `certifyMantleKernel` | {mat:func}`certifyShzKernel` |
| `geometry.flt2flt` | {mat:func}`geometry.splitSegment <splitSegment>` |
| `create_megathrust` | {mat:func}`createMegathrust` |

Error identifiers follow their functions, for example
`certifyShzKernel:notNegative`.

:::{warning}
`computePeriodicInitialCondition` existed before with the same arguments: it
was the legacy solver that replaced non-decaying modes by a fast decay. That
version is deleted and the name now belongs to the exact periodic solver
(formerly `computePeriodicInitialConditionBalanced`). A script that called the
old one still runs but now gets the exact solution, and raises an error if the
operator has growing or neutral modes instead of hiding them.
:::

## Deleted

| Function | Use instead |
| --- | --- |
| `computePeriodicInitialCondition` (legacy) | {mat:func}`computePeriodicInitialCondition` (exact) |
| `reconditionODEoperator` | the default certified kernels; for centre-sampled kernels, {mat:func}`negateGrowingModes` |
| `computeKLmodesMesh` | none (unused) |
| `computeSiteVelocitiesBem` | the displacement kernels of {mat:func}`computeAllDisplacementKernelsBem`, as in {doc}`examples/example1` |

## Other changes

- {mat:func}`createMegathrust` no longer leaves `megathrust2d.seg` in the
  current folder: it writes a temporary file and deletes it.
