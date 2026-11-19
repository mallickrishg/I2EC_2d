# Theory and methods

This page summarises the method and, just as important, the numerical
treatments that keep the discretised system stable. Read
[Caveats](#caveats-keeping-the-system-damped) before interpreting fine-scale
features of a solution. The full derivation is in Mallick & Sathiakumar (2026),
listed on the {doc}`contributors` page.

## Splitting the deformation

Deformation is split into two parts:

1. **Long-term component**: a quasi-static boundary value problem for the
   background flow (Newtonian corner flow, see
   {mat:func}`getStrainratesLongterm`) and steady fault slip at the plate rate
   $v^{\infty}$. At this timescale the elastic plates behave rigidly.
2. **Earthquake-cycle component**: a time-dependent perturbation about the
   long-term state, solved with an integral representation of the elastic
   response of the medium.

## Governing equations

The medium is discretised into fault patches on a boundary $\Gamma$ and
triangular cells (shear zones) in the mantle $\Omega$. The stress rate at any
point is

$$
\dot{\sigma}_i(\mathbf{x}) =
\int_{\Gamma} K_i(\mathbf{x};\xi)\,\bigl(v(\xi) - v^{\infty}\bigr)\,d\xi
+ \iint_{\Omega} K_{ji}(\mathbf{x};\alpha)\,
  \bigl(\dot{\epsilon}_{\mathrm{v}_j}(\alpha) -
        \dot{\epsilon}^{\infty}_{\mathrm{v}_j}(\alpha)\bigr)\,d\alpha ,
$$

where $K$ are elastostatic kernels that map fault slip and viscous strain to
stress. Only the deviatoric part of the mantle strain is evolved: the viscous
flow is incompressible, so $\dot\epsilon_{33} = -\dot\epsilon_{22}$, and the two
mantle unknowns per cell are $\dot\epsilon_{22}$ and $\dot\epsilon_{23}$.
After discretisation the kernels become the matrices returned by
{mat:func}`computeAllStressKernelsBem`:

| Kernel | Size | Maps |
| --- | --- | --- |
| `KK` | $N \times N$ | fault slip → fault shear traction |
| `KL` | $M \times N \times 2$ | fault slip → mantle deviatoric stress ($(\sigma_{22}-\sigma_{33})/2$, $\sigma_{23}$) |
| `LK` | $N \times M \times 2$ | mantle deviatoric strain → fault shear traction |
| `LL` | $M \times M \times 2 \times 2$ | mantle deviatoric strain → mantle deviatoric stress |

with $N$ fault patches and $M$ mantle cells. Units: MPa per metre of slip, or
MPa per unit strain.

## Rheology

**Mantle.** Each cell follows a power law,

$$
\dot{\epsilon}_{ij} = \alpha\, \tau^{\,n-1}\, \sigma_{ij},
\qquad \tau = \sqrt{\sigma_{22}^2 + \sigma_{23}^2},
$$

with the fluidity `shz.alpha` $= 1/\eta$ in the Newtonian case ($n = 1$) and
`shz.n` the stress exponent.

**Fault.** The two solvers treat the fault differently:

- *Numerical integration* uses rate-strengthening friction with parameter
  $(a-b)\sigma$ (`rcv.Asigma`, MPa), integrated in $\log v$:
  $\;d \log v/dt = \dot{\tau}/\bigl((a-b)\sigma\bigr)$.
- *Semi-analytical solution* replaces friction by a linear viscous fault,
  $dv/dt = \dot\tau/\eta'$, where `rcv.Asigma` now holds an effective
  viscosity per unit length $\eta' = \eta/L$ (MPa s/m). The example script
  uses $\eta = 5\times10^{17}$ Pa s spread over the creeping length.

Locked patches (`rcv.pinnedPosition`) have $v = 0$ between earthquakes and are
removed from the unknowns. Earthquakes are imposed: the slip deficit of the
locked zone, $v^{\infty} T$, is released every recurrence interval $T$, and its
stress change ($\Delta\tau = $ `KK`·slip on the fault, `KL`·slip in the
mantle) is applied as an instantaneous step.

## Solvers

### Numerical integration (general rheology)

{mat:func}`ode_viscoelastic` evaluates the right-hand side for slip rate
($\log v$) and mantle stress, and {mat:func}`viscoeq_sequence_simulation`
integrates it with MATLAB's adaptive Runge–Kutta 4(5) solver `ode45`
(absolute and relative tolerance $10^{-10}$), restarting after each imposed
event. {mat:func}`runImposedEarthquakeCycles` repeats the cycle `nreps` times
from the long-term state so that the last cycle is (close to) periodic; the
spin-up is brute force, and its convergence is not checked.

### Semi-analytical solution (linear case)

For $n = 1$ and a linear viscous fault, the state $y$ (creeping-fault slip rate
and mantle strain rates) obeys

$$
\frac{dy}{dt} = A\,y + c, \qquad A = \operatorname{diag}(1/\eta)\, K,
$$

where $K$ is assembled by {mat:func}`assembleStressKernel` (locked patches
removed), $\eta$ stacks the fault and cell viscosities, and $c$ is the forcing
from locking and long-term loading. With $A = V\,\mathrm{diag}(\lambda)\,V^{-1}$
each mode evolves on its own:

$$
z_k(t) = e^{\lambda_k t} z_k(0) + g_k(t)\, c_k,
\qquad g_k(t) = \frac{e^{\lambda_k t}-1}{\lambda_k},
$$

and the integrated state (slip, strain) uses
$h_k(t) = (g_k(t) - t)/\lambda_k$. {mat:func}`egh_modal` evaluates $e$, $g$
and $h$ exactly, with an 8-term Taylor series when $|\lambda t| < 10^{-2}$ so
that they stay accurate as $\lambda \to 0$. The solution at any time costs one
matrix–vector product ({mat:func}`analyticSolveIVP`).

**Periodic cycle.** For events with recurrence time $T$ that cause a jump $b$
in the state ($b = \Delta\sigma/\eta$), the steady cycle satisfies
$y(0) = y(T) + b$. In modal coordinates each mode gives
$z_k(0) = \bigl(g_k(T)\,c_k + b_k\bigr)/\bigl(1 - e^{\lambda_k T}\bigr)$
({mat:func}`computePeriodicInitialConditionBalanced`), so no spin-up is needed.
This is well defined only if every mode decays, $|e^{\lambda_k T}| < 1$, which
is what the treatments below guarantee.

(caveats-keeping-the-system-damped)=
## Caveats: keeping the system damped

Physically, $A$ describes a damped system: any perturbation relaxes back to the
long-term state, so every eigenvalue should have $\mathrm{Re}(\lambda) < 0$.
The discretised kernels do not quite honour this. Several treatments, applied
at different stages and in different combinations by the two solvers, make the
system damped again. They are regularisations, not physics.

### 1. How the kernels are computed

- **Half-space solutions plus a boundary correction.** Fault kernels use the
  plane-strain line-dislocation solutions of Davis & Martel (`LDstressHS`,
  `LDdispHS`); mantle kernels use the triangle strain-volume solutions of
  Barbot (2018). Both are for a half-space. The finite model box
  (`inputs/boundary2d.seg`) is imposed afterwards by a boundary-element solve
  for each source ({mat:func}`computeShzStressKernelsBem`,
  {mat:func}`computeFaultTractionKernelsBem`): prescribed displacement where
  `boundary.Vpl == 0` (Dirichlet) and zero traction where `boundary.Vpl == 1`
  (Neumann). There is no topography.
- **Boundary values must be zero.** The stress and traction kernel solves use
  `boundary.Vx`, `boundary.Vz` as the Dirichlet values, while the displacement
  kernel solves use zero. Kernels are responses to a unit source, so these
  values should be zero (as in the example, `bcx = bcz = 0`); a non-zero
  value would add the same constant to every column of the stress kernels.
- **Mantle stress by finite differences.** The stress kernels of a triangle are
  not evaluated from closed-form stresses: they are computed from the analytic
  *displacement* by finite differences with a step of $10^{-6}$ times the
  triangle's shortest edge (central in $x_2$, one-sided in $x_3$), at the cell
  centroids, with the source's own anelastic strain removed inside the source
  triangle ({mat:func}`computeStressPlaneStrainTriangleShearZoneFiniteDifference`).
  Where the closed-form displacement is not finite (on edges or vertices),
  tanh-sinh quadrature is used instead
  ({mat:func}`computeDisplacementPlaneStrainTriangleShearZoneTanhSinh`).
  This is the "regularised" stress kernel: it avoids the hypersingular
  closed-form evaluation, but self- and neighbour interactions remain the
  least accurate entries.
- **Deviatoric projection.** `KL`, `LK` and `LL` are combinations of the full
  $\sigma_{22}, \sigma_{23}, \sigma_{33}$ kernels with $\epsilon_{33} =
  -\epsilon_{22}$ imposed. The source code flags the index construction of
  the deviatoric `LL` kernel as still to be verified.

### 2. Spurious growing modes

The discretised `LL` kernel, and therefore $A$, has a few eigenvalues with
$\mathrm{Re}(\lambda) > 0$. They are associated with spatial scales comparable
to the cell size and appear on every mesh tried, including structured
meshes of uniform cells; remeshing does not remove them, so they have to be
handled numerically. Left alone, each would grow as $e^{\lambda t}$ and the
periodic condition above would be meaningless. The numerical (`ode45`) and
semi-analytical solvers handle them differently:

- **Change them in the kernel (numerical solves)** —
  {mat:func}`modifyGrowingModes`, also the `kernelmodify` option of
  {mat:func}`computeAllStressKernelsBem` (default `1`). Eigendecomposes the
  $2M \times 2M$ deviatoric `LL` kernel, changes the eigenvalues with
  $\mathrm{Re}(\lambda) > 0$ and rebuilds `LL`. Acts on the kernel alone,
  before viscosities are applied, so it works for any power-law exponent: for
  a reciprocal kernel the linearised mantle system has as many growing modes
  as `LL`, whatever $n$ and the stress state. Two choices:
  - `kernelmodify = 1` sets them to **zero**. The modes become neutral;
    nothing restores them, so they can accumulate over many earthquake
    cycles and slow the approach to a cycle-invariant state.
  - `kernelmodify = 2` sets them to $-\max|\lambda|$ (**damped**). They decay
    at the fastest rate the kernel has. Prefer this for long spin-ups and for
    power-law ($n = 3$) runs.

  *Used by* the numerical integration (`MAIN_i2ec_singleeqcycle.m`). Not for
  the semi-analytical solvers, which take the full kernel.
- **Negate them in the operator** — {mat:func}`negateGrowingModes`. Balances
  $A$, eigendecomposes it and replaces each $\lambda$ with
  $\mathrm{Re}(\lambda) > 0$ by $-\lambda$, keeping its eigenvector. Acts on
  the full coupled operator, after viscosities are applied. *Used by* the
  semi-analytical solver (`MAIN_analyticsolve_singleeqcycle.m`, through
  {mat:func}`analyticSolveIVP` and
  {mat:func}`computePeriodicInitialConditionBalanced`), with kernels computed
  with `kernelmodify = 0` (`kernels/evl_orig.mat`).
- **Replace them by a fast decay** — {mat:func}`reconditionODEoperator`. Sets
  eigenvalues with $\mathrm{Re}(\lambda) \ge 0$ to the real value
  $-0.1\ \mathrm{s^{-1}}$ (dropping their imaginary parts), so they decay
  within seconds; no balancing. *Used by* the legacy
  {mat:func}`computePeriodicInitialCondition` only.

Consequences to keep in mind:

- **Results at the cell scale are not reliable.** The modified modes are the
  mesh-scale ones, so the treatments are appropriate for smoothly varying
  strain-rate fields. Single-cell features and sharp cell-to-cell contrasts in
  viscosity or strain rate should not be interpreted.
- **The solvers are not interchangeable at the cell scale.** Zeroing an
  eigenvalue of `LL` leaves a neutral mode (it neither grows nor decays);
  negating it in $A$ makes it decay at the same rate it would have grown;
  replacing it makes it vanish within seconds. The two driver scripts and the
  legacy periodic solver can therefore give slightly different answers on the
  same mesh.
- **The set of modified modes depends on the viscosities.** `negateGrowingModes`
  works on $A = \operatorname{diag}(1/\eta) K$, whose eigenvalues change with
  $\eta$. Outputs vary continuously with viscosity, except where a
  complex-conjugate pair crosses $\mathrm{Re}(\lambda) = 0$: negation then
  flips the sign of its imaginary part, a jump of up to
  $2|\sin(\mathrm{Im}(\lambda)\,t)|$ of that mode's amplitude.
- **Forcing must use the treated operator.** Build every forcing term that
  contains the operator with the negated one, e.g.
  $c = \dot\epsilon_\text{lock} - A_\text{neg}\,y_\text{pl}$. With the raw $A$,
  the difference $(A_\text{neg} - A)\,y_\text{pl}$ acts as a spurious constant
  forcing and the plate-rate state is no longer steady.
- **Nearly neutral modes decay slowly.** A spurious eigenvalue close to zero
  stays close to zero after negation. Modes with $|\lambda| < 10^{-16}\
  \mathrm{s^{-1}}$ are treated as pure integrators and given zero periodic
  amplitude (the minimum-norm solution); this threshold is absolute, so the
  result is not continuous across it.

### 3. Other approximations

- **Interseismic steady state.** {mat:func}`computeInterseismicSteadystate`
  solves for the fully relaxed state with the fault locked. The relevant kernel
  block is singular, so the solution is regularised (Tikhonov, zeroth order,
  weight `lambda`; the notebook uses 10). The penalty pulls the mantle strain
  rates towards zero (not towards their long-term values) and leaves fault slip
  rates unpenalised, so the result depends on that weight. The function also
  assumes the locked patches are the first ones in `rcv`.
- **Corner-flow loading.** The long-term strain rates come from analytic
  isoviscous corner flow (Turcotte & Schubert, 2002), which is singular at the
  corner. The evaluation points are shifted (`ocean_translate_params`,
  `continental_translate_params`) to stay away from it, and the strain rate in
  each cell comes from a planar fit of the velocity at its three vertices. The
  shifts are tuning choices that should match the mesh.
- **Hinge fault.** Surface velocities and displacements include a vertical
  "hinge" fault (`inputs/hinge2d.seg`) that stands in for plate bending. It is
  kinematic: it is added to the surface predictions only and is not part of
  the dynamical system.
- **Imposed, characteristic earthquakes.** Events are prescribed (slip, timing),
  not computed; the periodic solution describes a perfectly repeating cycle.
- **Numerical safeguards in the linear solver.** $A$ is balanced (MATLAB
  `balance`) before its eigendecomposition to reduce round-off, and the series
  switch in {mat:func}`egh_modal` avoids cancellation in
  $(e^{\lambda t} - 1)/\lambda$ for small $|\lambda t|$.

### Practical guidance

- Use meshes fine enough that the features you interpret span several cells.
- After changing the mesh or the viscosity structure, check how many
  eigenvalues were modified and how large they were
  (`[~, lambda] = negateGrowingModes(A)` against `eig(A)`); a large number, or
  ones with large $|\lambda|$, suggests a mesh problem.
- Recompute the kernels after any geometry change; the files in `kernels/`
  match the example meshes only, and `evl_orig.mat` was computed with
  `kernelmodify = 0`. For a numerical solve, build the zeroed or damped kernel
  from it without recomputing: `evl = modifyGrowingModes(evl_orig, 2)`.
- `tests/test_negate_exact_solver.m` checks the semi-analytical solver
  against the matrix exponential, periodicity, continuity in viscosity and the
  steady state on the example mesh.
