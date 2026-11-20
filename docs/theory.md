# Theory and methods

This page summarises the model, the solvers and how the kernels are computed.
The full derivation is in Mallick & Sathiakumar (2026), listed on the
{doc}`contributors` page.

## Splitting the deformation

Deformation is split into two parts:

1. **Long-term component**: a quasi-static boundary value problem for the
   background flow (corner flow, Newtonian or power-law with $n = 3$, see
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

**Fault.** The fault law depends on the solver:

- *Numerical integration* uses, by default, rate-strengthening friction with
  parameter $(a-b)\sigma$ (`rcv.Asigma`, MPa), integrated in $\log v$:
  $\;d \log v/dt = \dot{\tau}/\bigl((a-b)\sigma\bigr)$
  (`rcv.faultLaw = 'ratestrengthening'`).
- *Semi-analytical solution* replaces friction by a linear viscous fault,
  $dv/dt = \dot\tau/\eta'$, where `rcv.Asigma` now holds an effective
  viscosity per unit length $\eta' = \eta/L$ (MPa s/m). The example script
  uses $\eta = 5\times10^{17}$ Pa s spread over the creeping length.
- The numerical integration can use the same linear viscous fault:
  set `rcv.faultLaw = 'linear'` (read by {mat:func}`isLinearFaultLaw`). The
  state is then $v$ instead of $\log v$, so $v$ may become negative. With
  $n = 1$ both solvers then integrate the same linear equations, which
  {doc}`examples/example4` uses to check one against the other.

Locked patches (`rcv.pinnedPosition`) have $v = 0$ between earthquakes and are
removed from the unknowns. Earthquakes are imposed: the slip deficit of the
locked zone, $v^{\infty} T$, is released every recurrence interval $T$, and its
stress change ($\Delta\tau = $ `KK`·slip on the fault, `KL`·slip in the
mantle) is applied as an instantaneous step.

## Solvers

### Numerical integration (general rheology)

{mat:func}`ode_viscoelastic` evaluates the right-hand side for slip rate
($\log v$, or $v$ for a linear viscous fault) and mantle stress, and {mat:func}`viscoeq_sequence_simulation`
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
from locking and long-term loading. $A$ is far from symmetric and badly
scaled (fault and mantle rows differ by their viscosities), so it is balanced
(MATLAB `balance`) before its eigendecomposition, which reduces round-off
({mat:func}`eigBalanced`). With $A = V\,\mathrm{diag}(\lambda)\,V^{-1}$
each mode evolves on its own:

$$
z_k(t) = e^{\lambda_k t} z_k(0) + g_k(t)\, c_k,
\qquad g_k(t) = \frac{e^{\lambda_k t}-1}{\lambda_k},
$$

and the integrated state (slip, strain) uses
$h_k(t) = (g_k(t) - t)/\lambda_k$. {mat:func}`egh_modal` evaluates $e$, $g$
and $h$ exactly, with an 8-term Taylor series when $|\lambda t| < 10^{-2}$ so
that they stay accurate for slow modes at short times, where the closed forms
lose digits to cancellation. The solution at any time costs one
matrix–vector product ({mat:func}`analyticSolveIVP`).

**Periodic cycle.** For events with recurrence time $T$ that cause a jump $b$
in the state ($b = \Delta\sigma/\eta$), the steady cycle satisfies
$y(0) = y(T) + b$. In modal coordinates each mode gives
$z_k(0) = \bigl(g_k(T)\,c_k + b_k\bigr)/\bigl(1 - e^{\lambda_k T}\bigr)$
({mat:func}`computePeriodicInitialConditionBalanced`), so no spin-up is needed.
This is well defined only if every mode decays, $|e^{\lambda_k T}| < 1$, which
the default mantle kernels guarantee (see [Mantle kernels](#mantle-kernels)).

(kernels)=
## Kernels

**Half-space solutions plus a boundary correction.** Fault kernels use the
plane-strain line-dislocation solutions of Davis & Martel (`LDstressHS`,
`LDdispHS`); mantle kernels use the triangle strain-volume solutions of
Barbot (2018). Both are for a half-space. The finite model box
(`examples/inputs/boundary2d.seg`) is imposed afterwards by a boundary-element solve
({mat:func}`computeShzStressKernelsBem`, {mat:func}`computeFaultTractionKernelsBem`):
zero displacement where `boundary.Vpl == 0` (Dirichlet) and zero traction where
`boundary.Vpl == 1` (Neumann). The boundary slip that cancels a mantle source
does not depend on where the stress is evaluated, so it is solved once for all
sources ({mat:func}`computeShzBoundaryCorrection`, one factorisation) and
reused for every receiver set. There is no topography.

**Homogeneous boundary problem.** All the BEM kernels (stress, traction and
displacement) have zero displacement on Dirichlet patches and zero traction on
Neumann patches. The prescribed boundary motion in the `.seg` file (`bcx`,
`bcz`) is not used: the background rate comes from the long-term solution.

**Analytic mantle stress.** The stress of a uniformly strained triangle is the
exact derivative of Barbot's (2018) closed-form displacement, generated with
the Symbolic Math Toolbox (`derivation/generateTriangleStress.m`), with the
source's own anelastic strain removed inside the triangle
({mat:func}`computeStressPlaneStrainTriangleShearZone`). Some terms, such as
$\arctan(N/D)/D$, cancel only between the two ends of an edge and lose
precision close to the lines through image edges and through horizontal or
vertical edges; there that edge's contribution is integrated numerically
(Gauss-Legendre or tanh-sinh quadrature). The stress agrees with an independent
quadrature to about $10^{-12}$, satisfies Betti reciprocity and the
traction-free surface to rounding, and is singular only on the triangle's own
edges, where cell-averaged kernels never evaluate it.

**Deviatoric projection.** `KL`, `LK` and `LL` are combinations of the full
$\sigma_{22}, \sigma_{23}, \sigma_{33}$ kernels with $\epsilon_{33} =
-\epsilon_{22}$ imposed. The source code flags the index construction of the
deviatoric `LL` kernel as still to be verified.

(mantle-kernels)=
### Mantle kernels: cell-averaged and certified

The continuous operator from a deviatoric eigenstrain to the stress it causes
is symmetric and negative semidefinite: stress always opposes the eigenstrain,
so viscous relaxation can only decay. Sampling the stress at each cell's
centroid (collocation) breaks both properties between neighbouring cells and
gives a few spurious growing modes (43 on the example mesh), however exactly
the stress is computed.

{mat:func}`computeAllStressKernelsBem` therefore averages the stress over each
receiving cell (3-point rule, option `cellaverage`). With piecewise-constant
strain this is a Galerkin discretisation, which keeps the symmetry and sign of
the continuous operator. The area-weighted `LL` is then symmetrised and checked
to be negative semidefinite ({mat:func}`certifyMantleKernel`, option
`kernelcertify`); the build fails if the check fails. A certified kernel has no
growing modes for any viscosity and any power-law exponent, so both solvers use
it as it is. The fault blocks (`KK`, `LK`) are not averaged.

**Centroid kernels** (`'cellaverage', 1, 'kernelcertify', false`) remain
available for comparison, with the earlier fixes for their growing modes:
`kernelmodify` 1 or 2 ({mat:func}`modifyGrowingModes`) sets them to zero or
damps them in `LL` before an `ode45` run, and {mat:func}`negateGrowingModes`
negates them in the operator $A$ for the semi-analytical solvers. These are
regularisations: results at the cell scale depend on the fix and should not be
interpreted, and with a negated operator every forcing term that contains it
(such as $c = \dot\epsilon_\text{lock} - A_\text{neg}\,y_\text{pl}$) must use
the negated operator too.

## Modelling assumptions and approximations

- **Corner-flow loading.** The long-term strain rates come from analytic
  corner flow (Turcotte & Schubert, 2002; $n = 1$ or $n = 3$), which is singular at the
  corner. The evaluation points are shifted (`ocean_translate_params`,
  `continental_translate_params`) to stay away from it, and the strain rate in
  each cell comes from a planar fit of the velocity at its three vertices. The
  shifts are tuning choices that should match the mesh.
- **Hinge fault.** Surface velocities and displacements include a vertical
  "hinge" fault (`examples/inputs/hinge2d.seg`) that stands in for plate bending. It is
  kinematic: it is added to the surface predictions only and is not part of
  the dynamical system.
- **Imposed, characteristic earthquakes.** Events are prescribed (slip, timing),
  not computed; the periodic solution describes a perfectly repeating cycle.

## Practical guidance

- Use meshes fine enough that the features you interpret span several cells
  ({doc}`examples/example3` shows the effect of mesh size).
- After building kernels for a new mesh, look at `evl.certificate`: the
  asymmetry removed and the largest eigenvalue of the symmetrised kernel. A
  failed certification is an error; it points to a mesh or kernel problem.
- Recompute the kernels after any geometry change. The files in `kernels/`
  match the example meshes only.
- `runTests` checks the Green's functions, kernels and solvers against
  physical identities; see `tests/README.md`.
