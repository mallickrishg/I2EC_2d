# Assumptions and valid range

What the model leaves out, and the range of mantle viscosities over which
each solver can be trusted.

## Assumptions

- **Imposed earthquakes.** Slip and timing are prescribed, not computed; the
  periodic solution is a perfectly repeating cycle.
- **Analytic corner flow.** The long-term loading is the corner-flow solution
  (Turcotte & Schubert, 2002), shifted away from its singular corner. The
  shifts (`ocean_translate_params`, `continental_translate_params`) are tuning
  choices that should match the mesh.
- **Kinematic hinge.** A vertical "hinge" fault (`examples/inputs/hinge2d.seg`)
  stands in for plate bending in the surface predictions only.

(valid-range)=
## Valid range of Maxwell times

How faithful a simulation is depends on the ratio of the mantle Maxwell time
$t_M = \eta/\mu$ to the recurrence time $T$, and on the fault, which sets the
fastest timescale. The two solvers have different limits.

**Semi-analytical solver.** It is accurate from very weak to very strong
mantles. The limit is reached only when the slowest mode relaxes by less than
$10^{-12}$ in one cycle ($|\lambda_\text{min}|\,T \le 10^{-12}$), when its
periodic state is no longer determined and the solver raises an error. On the
example 2 model with a uniform mantle:

| $t_M/T$ | $10^{-5}$ | $0.02$ (example 2) | $10^{2}$ | $10^{4}$ | $10^{6}$ | $10^{8}$ | $10^{9}$ | $10^{10}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| slowest $\lvert\lambda\rvert\,T$ | $3\times10^{2}$ | $0.14$ | $3\times10^{-5}$ | $3\times10^{-7}$ | $3\times10^{-9}$ | $3\times10^{-11}$ | $3\times10^{-12}$ | $3\times10^{-13}$ |
| relative difference, mantle state | $7\times10^{-11}$ | $4\times10^{-12}$ | $2\times10^{-10}$ | $5\times10^{-10}$ | $10^{-9}$ | $2\times10^{-9}$ | $6\times10^{-10}$ | error |

In the strong-mantle limit the mantle simply follows its long-term strain
rate, and the solver reproduces that.

**Numerical integration.** Accuracy is set by the `ode45` tolerances, but each
cycle only brings the solution closer to the periodic state by the factor
$e^{-|\lambda_\text{min}|\,T}$. Reaching $10^{-6}$ takes about
$14/(|\lambda_\text{min}|\,T)$ cycles: a few for $t_M/T \lesssim 10^{-3}$, about
100 for example 2, about 5000 for $t_M/T = 1$ and about $5\times10^{5}$ for
$t_M/T = 10^{2}$. For linear rheology beyond $t_M/T \sim 0.1$, use the
semi-analytical solver. For power-law rheology the effective Maxwell time
changes through the cycle, so check convergence directly, as example 1 does.

Both limits shift with the fault: a stiffer or weaker fault changes the
fastest timescale and so the numbers above.

:::{admonition} In the code
:class: in-code

the table compares {mat:func}`computePeriodicInitialCondition` with the more
accurate of two independent eigendecomposition solves, balanced for weak
mantles and unbalanced for strong ones (example 2 fault, locked above 30 km,
$T = 500$ yr). Before
2026-11-22 the periodic state came from a single eigendecomposition: balanced
(until 2026-11-21) it lost digits for strong mantles ($10^{-8}$ at
$t_M/T = 10^{4}$, $10^{-6}$ at $10^{6}$), unbalanced it lost them for the
example model ($3\times10^{-9}$ in the mantle state). Modes that relax by less than $10^{-12}$ per cycle
($|\mathrm{Re}\,\lambda|\,T \le 10^{-12}$) raise `I2EC:neutralModes`, since
their periodic state is not determined. Before 2026-11-21 the solver instead
set modes with $|\lambda| < 10^{-16}$ s$^{-1}$ to zero, which for
$t_M/T \gtrsim 10^{4}$ (with this fault) dropped real mantle modes and left the
periodic mantle state 58 % to 100 % wrong: recompute strong-mantle results made
with earlier versions. `ode45` is explicit, so its steps per cycle also grow
with the fastest rate, $\max|\lambda|\,T$; with this fault the fault sets it
down to $t_M/T \approx 10^{-5}$.
:::

## Practical guidance

- Check that the Maxwell time of your model is in the range the solver
  handles ({ref}`valid-range`).
- Use meshes fine enough that the features you interpret span several cells;
  {doc}`../examples/example3` shows the effect of mesh size.
- Recompute the kernels after any geometry change; the files in
  `examples/kernels/` match the example meshes only.
- After building kernels for a new mesh, look at `evl.certificate`. A failed
  certification is an error that points to a mesh or kernel problem.
- `runTests` checks the Green's functions, kernels and solvers against
  physical identities; see `tests/README.md`.
