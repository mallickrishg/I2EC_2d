# Theory and methods

I2EC models the earthquake cycle as small departures from steady plate motion.
The fault and the mantle are discretised into patches and triangular cells,
elastic kernels link them, and the result is a system of ordinary differential
equations that is either integrated in time or, when it is linear, solved
exactly. The full derivation is in Mallick & Sathiakumar (2026), listed on the
{doc}`contributors` page.

**In five steps:**

1. [Long-term loading](#long-term-loading): steady plate motion and mantle
   corner flow.
2. [Governing equations](#governing-equations): fault slip and mantle flow
   change the stress through elastic kernels.
3. [Rheology](#rheology): the stress sets the fault slip rate and the mantle
   strain rate.
4. [Solvers](#solvers): integrate in time, or solve the linear case exactly.
5. [Kernels](#kernels): how the elastic kernels are computed.

## Long-term loading

Deformation is split into a **long-term** part and an **earthquake-cycle**
part. Over geological time the plates move steadily: the fault slips at the
plate rate $v^{\infty}$ and the mantle flows as corner flow (Newtonian, or
power-law with $n = 3$). The cycle is a time-dependent departure from this
state, driven by locking and imposed earthquakes.

:::{admonition} In the code
:class: in-code

{mat:func}`computeLongTermStrainRates` gives the corner-flow
strain rates in each cell, {mat:func}`computeLongTermVelocity` the velocity
field ({doc}`examples/index` plots it).
:::

## Governing equations

The fault is a set of patches on a boundary $\Gamma$ and the mantle a set of
triangular cells $\Omega$. Slip and viscous strain that depart from the
long-term state change the stress everywhere:

$$
\dot{\sigma}_i(\mathbf{x}) =
\int_{\Gamma} K_i(\mathbf{x};\xi)\,\bigl(v(\xi) - v^{\infty}\bigr)\,d\xi
+ \iint_{\Omega} K_{ji}(\mathbf{x};\alpha)\,
  \bigl(\dot{\epsilon}_{\mathrm{v}_j}(\alpha) -
        \dot{\epsilon}^{\infty}_{\mathrm{v}_j}(\alpha)\bigr)\,d\alpha .
$$

The kernels $K$ are elastic Green's functions. Mantle flow is incompressible,
so each cell has two unknowns, the deviatoric strain rates
$\dot\epsilon_{22}$ and $\dot\epsilon_{23}$. Discretised, the kernels are four
matrices ($N$ fault patches, $M$ mantle cells):

| Kernel | Maps | Size |
| --- | --- | --- |
| `KK` | fault slip → fault shear traction | $N \times N$ |
| `KL` | fault slip → mantle stress | $M \times N \times 2$ |
| `LK` | mantle strain → fault shear traction | $N \times M \times 2$ |
| `LL` | mantle strain → mantle stress | $M \times M \times 2 \times 2$ |

:::{admonition} In the code
:class: in-code

{mat:func}`computeAllStressKernelsBem` returns all four, in MPa
per metre of slip or per unit strain. Mantle stresses are the deviatoric pair
$(\sigma_{22}-\sigma_{33})/2$ and $\sigma_{23}$.
:::

## Rheology

**Mantle.** Each cell flows as a power law in the deviatoric stress $\tau$:

$$
\dot{\epsilon}_{ij} = \alpha\, \tau^{\,n-1}\, \sigma_{ij},
\qquad \tau = \sqrt{\sigma_{22}^2 + \sigma_{23}^2}.
$$

With $n = 1$ the mantle is Newtonian and $\alpha = 1/\eta$.

**Transient creep (Burgers, Andrade).** Each cell may also have Kelvin
elements, a spring $\mu_K$ and a dashpot $\eta_K$ in parallel, in series with
the power-law (Maxwell) dashpot:

$$
\sigma_{ij} = \mu_{K,k}\, \epsilon^{K}_{ij,k} + \eta_{K,k}\, \dot\epsilon^{K}_{ij,k},
\qquad
\dot\epsilon_{ij} = \alpha\,\tau^{\,n-1}\sigma_{ij} + \sum_k \dot\epsilon^{K}_{ij,k}.
$$

One element per cell is a Burgers body. Andrade creep, $\epsilon(t) = \beta
\sigma t^{1/3}$ under constant stress, is approximated by a Prony series:
Kelvin elements with log-spaced relaxation times whose compliances are fitted
to $t^{1/3}$ (about 1% with one term per decade). The kernels see only the
total strain rate, so the long-term flow is that of the Maxwell element, and
over a periodic cycle the Kelvin strain rates average to zero.

**Fault.** The creeping fault follows one of two laws:

| Law | Equation | Used by |
| --- | --- | --- |
| rate-strengthening friction | $\dfrac{d \log v}{dt} = \dfrac{\dot{\tau}}{(a-b)\sigma}$ | numerical solver (default) |
| linear viscous | $\dfrac{dv}{dt} = \dfrac{\dot\tau}{\eta'}$ | both solvers |

**Locking and earthquakes.** Locked patches do not slip between earthquakes.
Every recurrence interval $T$ an earthquake releases their slip deficit,
$v^{\infty} T$, and its stress change is applied as an instantaneous step.

:::{admonition} In the code
:class: in-code

the mantle takes `shz.alpha` and `shz.n`, and Kelvin elements in
`shz.etaK` (MPa s; `Inf` for none) and `shz.muK` (MPa), $M \times p$ each
(default $M \times 0$: Maxwell). {mat:func}`computeAndradePronySeries`
returns them for Andrade creep; {doc}`examples/example5` compares the three
rheologies. The fault takes
`rcv.faultLaw` (`'ratestrengthening'` or `'linear'`, read by
{mat:func}`isLinearFaultLaw`) and `rcv.Asigma`, which holds $(a-b)\sigma$ (MPa)
or $\eta'$, a viscosity per unit length (MPa s/m). Locked patches are
`rcv.pinnedPosition`.
:::

## Solvers

| | Numerical | Semi-analytical |
| --- | --- | --- |
| Rheology | any (power law, friction; Burgers, Andrade) | linear: $n = 1$ (Maxwell, Burgers, Andrade), linear viscous fault |
| Method | time integration with `ode45` | eigendecomposition, exact in time |
| Periodic cycle | spin-up over many cycles | solved directly |
| Example | {doc}`examples/example1` | {doc}`examples/example2` |

{doc}`examples/example4` runs both on the same linear model; they agree to
about $10^{-9}$.

### Numerical integration

The equations are integrated with MATLAB's adaptive Runge–Kutta solver
`ode45` (tolerance $10^{-10}$), restarting after each earthquake. Repeating the
cycle from the long-term state spins the model up until it is periodic.

:::{admonition} In the code
:class: in-code

{mat:func}`runImposedEarthquakeCycles` repeats the cycle;
{mat:func}`integrateImposedCycles` integrates one sequence of events and
{mat:func}`odeViscoelastic` is the right-hand side. Spin-up convergence is
not checked: compare the start and end of the last cycle, as example 1 does.
Kelvin elements add their spring stresses $\mu_K\epsilon^K$ to the state;
short relaxation times make the system stiff, so both drivers take an
optional solver such as `@ode15s`.
:::

### Semi-analytical solution

When everything is linear, the state $y$ (slip rate on the creeping patches
and mantle strain rates) obeys

$$
\frac{dy}{dt} = A\,y + c, \qquad A = \operatorname{diag}(1/\eta)\, K,
$$

with $c$ the forcing from locking and long-term loading. In the eigenvectors of
$A$ each mode evolves on its own,

$$
z_k(t) = e^{\lambda_k t} z_k(0) + \frac{e^{\lambda_k t}-1}{\lambda_k}\, c_k ,
$$

so the solution at any time costs one matrix–vector product. If each
earthquake makes the state jump by $b$, the periodic cycle,
$y(0) = y(T) + b$, has the closed form

$$
z_k(0) = \frac{g_k(T)\,c_k + b_k}{1 - e^{\lambda_k T}},
\qquad g_k(T) = \frac{e^{\lambda_k T}-1}{\lambda_k},
$$

which needs every mode to decay, $\mathrm{Re}(\lambda_k) < 0$. The default
kernels guarantee this ([mantle kernels](#mantle-kernels)).

:::{admonition} In the code
:class: in-code

{mat:func}`assembleStressKernel` builds $K$ without the
locked patches, and {mat:func}`assembleLinearSystem` builds $A$, $c$ and $b$.
With Kelvin elements it appends their strain rates to $y$: each row is
divided by the viscosity of its own element, $K$ acts on the total rates, and
Kelvin rows decay at $\mu_K/\eta_K$, so $A$, $c$ and $b$ go unchanged into
the periodic solver; {mat:func}`computePeriodicInitialCondition` gives the
periodic state and {mat:func}`solveAnalyticIvp` evaluates it at any time.
The fault rows of $A$ are about $10^6$ times larger than the mantle rows, so
no single eigendecomposition is accurate for every mantle: with balancing the
slow modes of a strong mantle lose digits, without it the mantle part of a
weak or moderate one does. {mat:func}`computePeriodicInitialCondition`
therefore solves the periodic condition directly, with one matrix exponential
of $\begin{bmatrix} A & I \\ 0 & 0 \end{bmatrix} T$, and uses the
eigendecomposition ({mat:func}`eigOperator`) only to check for growing and
neutral modes. {mat:func}`solveAnalyticIvp` eigendecomposes $A$ with and
without balancing and keeps the one that propagates the state more accurately,
checked against the matrix exponential. {mat:func}`computeModalIntegrals`
switches to a Taylor series for small $|\lambda t|$, where the closed forms
lose digits.
:::

(kernels)=
## Kernels

The kernels are **half-space** solutions with a **boundary correction**:

- **Fault patches** use the plane-strain line-dislocation solutions of Davis &
  Martel.
- **Mantle cells** use the analytic stress of a uniformly strained triangle:
  the exact derivative of Barbot's (2018) closed-form displacement. It agrees
  with an independent quadrature to about $10^{-12}$ and satisfies Betti
  reciprocity to rounding.
- **The model box** is imposed by a boundary-element solve: zero displacement
  on its fixed sides and zero traction on its free sides. The background rate
  comes from the long-term solution, not from boundary motion.

:::{admonition} In the code
:class: in-code

{mat:func}`computeStressPlaneStrainTriangleShearZone` is the
triangle stress, generated with the Symbolic Math Toolbox by
`derivation/generateTriangleStress.m`; near a few lines through its edges,
where the closed form loses precision, that edge's contribution is integrated
numerically. {mat:func}`computeShzBoundaryCorrection` solves the boundary
problem once for all mantle sources. The fixed and free sides are set by
`boundary.Vpl` (0 and 1) in `examples/inputs/boundary2d.seg`.
:::

(mantle-kernels)=
### Mantle kernels: cell-averaged and certified

Physically, the stress caused by a viscous strain always opposes it, so
viscous relaxation can only decay. Sampling the stress at each cell's centre
breaks this between neighbouring cells and creates a few spurious growing
modes (43 on the example mesh), however exactly the stress is computed.

I2EC therefore **averages the stress over each cell**, a Galerkin
discretisation that keeps the right symmetry and sign, and **certifies** the
result: the mantle–mantle kernel is checked to have no growing mode and no
zero-energy (neutral) mode. A certified kernel is stable for any viscosity and
any power-law exponent, in both solvers.

:::{admonition} In the code
:class: in-code

this is the default of {mat:func}`computeAllStressKernelsBem`
(options `cellaverage` and `kernelcertify`, check by
{mat:func}`certifyShzKernel`, recorded in `evl.certificate`). Centre-sampled
kernels (`'cellaverage', 1, 'kernelcertify', false`) remain for comparison.
Their growing modes need a fix: `kernelmodify` ({mat:func}`modifyGrowingModes`)
for `ode45`, or {mat:func}`negateGrowingModes` on $A$, with every forcing term
built from the negated operator. Their results at the cell scale should not be
interpreted.
:::

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
$T = 500$ yr). For $t_M/T \gtrsim 10^{6}$ MATLAB may warn that a matrix is
close to singular; the result still agrees to about $10^{-9}$. Before
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
  {doc}`examples/example3` shows the effect of mesh size.
- Recompute the kernels after any geometry change; the files in
  `examples/kernels/` match the example meshes only.
- After building kernels for a new mesh, look at `evl.certificate`. A failed
  certification is an error that points to a mesh or kernel problem.
- `runTests` checks the Green's functions, kernels and solvers against
  physical identities; see `tests/README.md`.
