# Example 6: viscosity inversion with adjoint gradients

A synthetic inversion: recover a smooth mantle viscosity $\eta(x)$ of the
example model from 20 years of surface displacement after an earthquake. The
functions in `inverse/` give the gradient of the misfit with respect to
$\log\eta$ in every mantle cell for about the cost of one forward solve, and
the Gauss–Newton Jacobian of a low-rank parameterisation from a few tangent
solves. Run it with

```matlab
run examples/example6_adjointInversion
```

It takes under a minute; each Gauss–Newton iteration (forward, 6 tangents)
takes about 0.6 s.

## 1. The forward model

The model is linear (Maxwell here; Burgers works the same way), so it is
solved without time stepping:

1. **Spin-up:** the periodic cycle of {doc}`example2`, one earthquake every
   500 years, solved directly (no spin-up cycles).
2. **Observed window:** that cycle's next earthquake at $t = 0$, then 40
   epochs over 20 years of horizontal and vertical displacement at 21 surface
   stations.

{mat:func}`solveLinearSequence` solves both, given the linear system of
{mat:func}`assembleLinearSystem` and the time line:

```matlab
sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change);
seq = struct('T', T, 'b0', sys.b, 'age', T, ...         % spin-up cycle
             'tEvents', 0, 'bEvents', sys.b, ...          % observed events
             'tObs', tObs, 'tRef', tObs(1));              % epochs, reference
fwd = solveLinearSequence(sys, seq);
```

`fwd.U` holds the slip and strain accumulated from `tRef` to each epoch.

## 2. Three layers

| Layer | Functions | Works on |
| --- | --- | --- |
| core | {mat:func}`computeAdjointGradient`, {mat:func}`computeTangentSequence` | state rows: $\partial J/\partial\log\eta$ per row, given $\partial J/\partial U$ |
| parameters | {mat:func}`extractViscosityGradient`, {mat:func}`assembleRowPerturbation` | fault patches and mantle cells (Maxwell: $\eta_M$; Burgers: $\eta_M$, $\tau_K$, $\mu_K/G$) |
| observations | {mat:func}`assembleStationKernels`, {mat:func}`computeDisplacementMisfit` | stations, epochs, weights, a free offset per station |

The fourth layer, the parameterisation, belongs to the caller. Here it is
$\log\eta = \log(10^{19}\,\mathrm{Pa\,s}) + \Phi\xi$, with $\Phi$ the first
6 Karhunen–Loève modes of a Matérn ($\nu = 3/2$, 60 km) covariance on the
mesh, so the gradient is $\Phi^\mathsf{T}\nabla$ and each Jacobian column is a
tangent along one mode:

```matlab
obs.H = assembleStationKernels(devl, sys, 'Stations', sta);
obs.d = d; obs.W = 1/sigma^2; obs.fixed = fixed; obs.offsets = true;
[J, res] = computeDisplacementMisfit(obs, fwd);
g  = computeAdjointGradient(fwd, res.S);        % per state row
gp = extractViscosityGradient(g, sys);          % per cell: gp.etaM
gradXi = Phi' * gp.etaM;

dTheta = assembleRowPerturbation(sys, struct('etaM', Phi));
dU = computeTangentSequence(fwd, dTheta);       % 6 tangents
[J, res] = computeDisplacementMisfit(obs, fwd, dU);   % res.Jw: Jacobian
```

`obs.fixed` is the part of the displacement that does not depend on
viscosity: here the long-term motion and the kinematic hinge. The script
checks the adjoint gradient against the tangent Jacobian
($J_w^\mathsf{T} r_w$):

```text
adjoint vs tangent gradient: relative difference 2.9e-16
```

## 3. Gauss–Newton

From a uniform $10^{19}$ Pa s mantle, Levenberg–Marquardt steps on $\xi$
with a unit Gaussian prior reach the noise level (5 mm) in three iterations:

```text
iteration 0: chi^2/N = 308.753, |xi - xi_true| = 1.715
iteration 1: chi^2/N =  15.575, |xi - xi_true| = 0.484
iteration 2: chi^2/N =   1.193, |xi - xi_true| = 0.059
iteration 3: chi^2/N =   1.032, |xi - xi_true| = 0.018
xi true:      [-1.2 0.8 0.6 -0.5 0.4 -0.3]
xi recovered: [-1.2 0.8 0.59 -0.5 0.39 -0.3]
```

```{figure} ../_static/examples/example6_adjointInversion_1.png
:alt: True and recovered mantle viscosity on the mesh
:width: 100%

True and recovered $\log_{10}(\eta/10^{19}\,\mathrm{Pa\,s})$. Triangles mark
the stations.
```

```{figure} ../_static/examples/example6_adjointInversion_2.png
:alt: Misfit per iteration and displacement time series at three stations
:width: 100%

Left: data misfit per iteration. Right: horizontal displacement at three
stations: data (dots), the starting model (dashed) and the final model
(solid).
```

## How the gradient is computed

Every row of the operator, the forcing and the coseismic jump is divided by
the viscosity of its element, so the derivative with respect to
$\log\eta_r$ is $-\int \lambda_r\,\mathrm{d}y_r$: the adjoint $\lambda$
against the change of the state, jumps included. The periodic spin-up adds
the adjoint of the periodicity condition, which is diagonal in modal form.
Both are closed-form in the modes of $A$, which the forward has already
computed, so the gradient costs one more $n \times n$ matrix product. The
derivation, the tests and the cost are in `docs/dev/adjoint.md`.
