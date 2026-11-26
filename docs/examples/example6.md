# Example 6: Burgers mantle inversion with adjoint gradients

A synthetic inversion: recover a Burgers mantle whose three properties all
vary in space (steady-state viscosity $\eta_M$, Kelvin viscosity $\eta_K$ and
Kelvin modulus $\mu_K$) from 20 years of surface displacement after an
earthquake. The functions in `inverse/` give the gradient of the misfit with
respect to every cell's parameters for about the cost of one forward solve,
and the Gauss–Newton Jacobian of a low-rank parameterisation from one tangent
solve per coefficient. Run it with

```matlab
run examples/example6_adjointInversion
```

It takes about 1.5 minutes. Each Gauss–Newton iteration (forward plus 12
tangents, 2103 states) takes about 5 s.

## 1. The forward model

The model is linear, so it is solved without time stepping:

1. **Spin-up:** the periodic cycle of {doc}`example2`, one earthquake every
   500 years, solved directly (no spin-up cycles).
2. **Observed window:** that cycle's next earthquake at $t = 0$, then 20 years
   of horizontal and vertical displacement at 21 surface stations every 3
   months (5 mm noise, plus an unknown offset per station and component).

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

## 2. Parameters: logs of ratios

Each cell has three Burgers parameters. They are not inverted as raw
$(\eta_M, \eta_K, \mu_K)$: the data see the Kelvin element mainly through its
relaxation time and its strength, so $\eta_K$ and $\mu_K$ alone trade off
almost perfectly. The parameters are

$$
\log\eta_M, \qquad \log\tau_K = \log(\eta_K/\mu_K), \qquad \log(\mu_K/G).
$$

Each is a reference value ($10^{19}$ Pa s; $\eta_K = 10^{18}$ Pa s and
$\mu_K = G$, so $\tau_K \approx 1$ yr) plus a few Karhunen–Loève modes of a
Matérn ($\nu = 3/2$, 60 km) covariance on the mesh: 5, 4 and 3 modes, with
prior standard deviations 1, 1 and 0.3 in natural log. These modes are the
caller's layer, built in the script.

To **hold $\mu_K$ fixed** and invert only $\eta_M$ and $\eta_K$, leave the
$\mu_K/G$ field out of the directions. With $\mu_K$ fixed, a change of
$\log\tau_K$ is the same change of $\log\eta_K$, so the ratio
parameterisation still holds. Case B does this.

## 3. Three layers in `inverse/`

| Layer | Functions | Works on |
| --- | --- | --- |
| core | {mat:func}`computeAdjointGradient`, {mat:func}`computeTangentSequence` | state rows: $\partial J/\partial\log\eta$ (and $\log\mu_K$) per row, given $\partial J/\partial U$ |
| parameters | {mat:func}`extractViscosityGradient`, {mat:func}`assembleRowPerturbation` | fault patches and mantle cells: $\eta_M$, $\tau_K$, $\mu_K/G$ |
| observations | {mat:func}`assembleStationKernels`, {mat:func}`computeDisplacementMisfit` | stations, epochs, weights, a free offset per station |

With `PhiM`, `PhiT` and `PhiU` the mode matrices of the three fields, the
gradient and the Jacobian are

```matlab
obs.H = assembleStationKernels(devl, sys, 'Stations', sta);
obs.d = d; obs.W = 1/sigma^2; obs.fixed = fixed; obs.offsets = true;

% gradient: adjoint per row -> per cell -> per coefficient
[J, res] = computeDisplacementMisfit(obs, fwd);
gp = extractViscosityGradient(computeAdjointGradient(fwd, res.S), sys);
gradXi = [PhiM'*gp.etaM; PhiT'*gp.tauK; PhiU'*gp.muKrel];

% Jacobian: one tangent per coefficient (dp fields: [M x r] or [M x 1 x r])
dp.etaM = ...; dp.tauK = ...; dp.muKrel = ...;   % leave muKrel out to fix mu_K
[dTheta, dPhi] = assembleRowPerturbation(sys, dp);
dU = computeTangentSequence(fwd, dTheta, dPhi);
[J, res] = computeDisplacementMisfit(obs, fwd, dU);   % res.Jw, res.rw
```

`obs.fixed` is the part of the displacement that does not depend on the
mantle: here the long-term motion and the kinematic hinge. The script checks
the adjoint gradient against the tangent Jacobian ($J_w^\mathsf{T} r_w$) in
both cases:

```text
adjoint vs tangent gradient: relative difference 1.1e-15
```

## 4. Gauss–Newton

From the reference mantle, Levenberg–Marquardt steps with a unit Gaussian
prior on the coefficients reach the noise level in three iterations, with
$\mu_K$ free (case A, 12 coefficients) and fixed at its known field (case B,
9 coefficients):

```text
case A: eta_M, tauK, mu_K/G free (12 coefficients)
final: chi^2/N = 1.023, |xi - xi_true| = 0.084
  xi true:      [-1 0.7 0.5 -0.4 0.3 0.8 -0.6 0.5 -0.3 1 -0.7 0.5]
  xi recovered: [-0.98 0.7 0.5 -0.41 0.29 0.83 -0.6 0.51 -0.31 0.99 -0.77 0.47]

case B: mu_K fixed, eta_M and tauK = eta_K/mu_K free (9 coefficients)
final: chi^2/N = 1.024, |xi - xi_true| = 0.046
```

```{figure} ../_static/examples/example6_adjointInversion_1.png
:alt: True and recovered Burgers parameters on the mesh, for both cases
:width: 100%

The three Burgers fields relative to their reference values: true (left),
recovered with all three free (middle) and with $\mu_K$ fixed (right).
```

```{figure} ../_static/examples/example6_adjointInversion_2.png
:alt: Misfit per iteration and displacement time series at three stations
:width: 100%

Left: data misfit per iteration for both cases. Right: horizontal
displacement at three stations in case A: data (dots), the reference mantle
(dashed) and the final model (solid). The fast early transient is the Kelvin
element's.
```

## How the gradient is computed

Every row of the operator, the forcing and the coseismic jump is divided by
the viscosity of its element (fault patch, Maxwell or Kelvin dashpot), so the
derivative with respect to $\log\eta_r$ is $-\int \lambda_r\,\mathrm{d}y_r$:
the adjoint $\lambda$ against the change of the state, jumps included. The
Kelvin spring enters only the diagonal, so the derivative with respect to
$\log\mu_K$ is $-(\mu_K/\eta_K)\int\lambda_r y_r\,\mathrm{d}t$. The periodic
spin-up adds the adjoint of the periodicity condition, which is diagonal in
modal form. All of these are closed-form in the modes of $A$, which the
forward has already computed, so the gradient costs one more $n \times n$
matrix product. The derivation, the tests and the cost are in
`docs/dev/adjoint.md`.
