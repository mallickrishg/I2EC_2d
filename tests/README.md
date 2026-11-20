# I2EC tests

A MATLAB `matlab.unittest` suite that checks what the code gets right and records what it gets wrong. Run it before and after every clean-up step in [`tasks.md`](../tasks.md).

## Running

From the repository root, in MATLAB:

```matlab
addpath tests
runTests            % fast tier, ~1-2 min (first run also starts a parallel pool)
runTests('all')     % + slow tier: rebuilds the example kernels, runs the MAIN scripts
```

Headless:

```bash
matlab -batch "addpath tests; runTests"
```

`runTests` prints a summary and writes `tests/results/summary.md` and `tests/results/junit.xml` (git-ignored). It throws an error, so `matlab -batch` exits non-zero, **only on unexpected failures**.

To run one class or one test:

```matlab
addpath tests/helpers tests/unit
runtests('KernelsTest')
runtests('KernelsTest/faultMantleKernelsAreReciprocal')
```

## How to read the results

| Group | Meaning |
|---|---|
| passed | the behaviour is correct |
| unexpected failure | a test **without** the `KnownIssue` tag failed: something that used to work is broken |
| known issue | a test tagged `KnownIssue` failed, as expected: it checks the *correct* behaviour of code that is currently wrong (listed in `tasks.md`) |
| known issue now passing | the code was fixed: remove the `KnownIssue` tag so the test guards the fix |

Tests tagged `Slow` only run with `runTests('all')`.

## Layout

```
tests/
  runTests.m                  runner and summary
  helpers/                    path fixture, synthetic and example models, small utilities
  unit/
    GreensFunctionsTest.m     line dislocations, triangle strain volumes (Barbot 2018)
    GeometryTest.m            +geometry classes, mesh discretisation, half-space kernels
    KernelsTest.m             BEM stress/displacement kernels, Betti reciprocity, cached kernels
    LongTermFlowTest.m        corner flow (n = 1, 3), long-term strain rates and velocity
    SolversTest.m             egh_modal, negate + modal solvers, periodic state, helpers
    TimeIntegrationTest.m     ode_viscoelastic / ode45 path vs the semi-analytical solver
    UtilitiesTest.m           plotting and small helpers
    LegacyCodeTest.m          duplicate copies, dead and broken files, name checks
  system/
    ExampleScriptsTest.m      end-to-end runs of examples/ (Slow)
  test_negate_exact_solver.m  stand-alone verification of the negate + exact solver
                              (also run by ExampleScriptsTest)
```

## What the checks rely on

Most tests check physical identities that hold whatever the sign convention, so a passing test means more than "the code runs":

- **Free surface**: stresses from every Green's function have zero traction at z = 0.
- **Dislocation jump**: displacement jumps by exactly the prescribed slip across a fault.
- **Hooke consistency**: stress kernels equal Hooke's law applied to the gradient of the displacement kernels.
- **Betti reciprocity**: weighted by patch length and triangle area, the KK, KL/LK and LL kernel blocks are symmetric. A sign or component-ordering error in any block breaks this.
- **Corner flow**: velocities satisfy their boundary conditions and are divergence-free.
- **Two solution routes agree**: with a locked fault and a Newtonian mantle the ODE path (`ode45`) is linear and must match `analyticSolveIVP`, through a coseismic step.
- **References**: the semi-analytical solvers are checked against `expm` of the augmented system.

Fast tests use a 10-patch / 24-triangle synthetic model (`helpers/makeSmallModel.m`), so the full BEM pipeline runs in about a second. Tests on the example mesh use `examples/inputs/` and the cached `kernels/*.h5`.

## Adding a test

Put it in the class for its layer and use `I2ECPathFixture` as a shared fixture, as the existing classes do. Name the test after the behaviour it checks, for example `periodicStateClosesTheCycle`. If the test documents a bug, put it in a `methods (Test, TestTags = {'KnownIssue'})` block, say in a comment what is wrong, and add the item to `tasks.md`.
