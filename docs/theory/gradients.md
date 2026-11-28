# Adjoint gradients

For linear models (Maxwell or Burgers mantle, linear viscous or Burgers
fault), I2EC computes the gradient of a data misfit with respect to the
viscosity of every mantle cell and fault patch, for about the cost of one
forward solve. This is what an inversion of surface displacement for mantle
viscosity needs. Power-law rheology has no adjoint.

## Misfit

The forward model is the semi-analytical solution of
{doc}`solvers`: a periodic spin-up cycle, then an observed window of
earthquakes and observation epochs. Surface displacement at the stations
is a linear function of the time-integrated state,

$$
u(t_j) = H \int_{t_\mathrm{ref}}^{t_j} y\,dt + u_\mathrm{fixed}(t_j) + B\,o,
\qquad
J = \tfrac12 \sum_j r_j^\mathsf{T} W r_j,\quad r_j = u(t_j) - d_j ,
$$

with $H$ the displacement kernels, $u_\mathrm{fixed}$ the part that does not
depend on viscosity (long-term motion, the kinematic hinge) and $o$ a free
offset per station and component, eliminated from the misfit.

## Adjoint

Every row of $A$, of the forcing $c$ and of the coseismic jump $b$ is divided
by the viscosity of its own element (fault patch, Maxwell or Kelvin dashpot).
For $\theta_r = \log\eta_r$ the derivative is therefore

$$
\frac{\partial J}{\partial \theta_r} = -\int \lambda_r\, dy_r ,
\qquad
-\dot\lambda = A^\mathsf{T}\lambda + s(t),
$$

the adjoint $\lambda$ against the change of the state, earthquake jumps
included. The adjoint runs backward in time, driven by the weighted
residuals $s(t)$, which are constant between epochs. The periodic spin-up
adds the adjoint of the periodicity condition, which is diagonal in the
eigenmodes. Every term is a closed-form integral of exponentials in the
eigenmodes of $A$ that the forward solve has already computed, so the
gradient for all rows costs one more $n \times n$ matrix product.

The **tangent** (forward sensitivity) along a direction $\delta\theta$,
$\delta\dot y = A\,\delta y - \delta\theta \circ \dot y$, uses the same
modal integrals. One tangent per parameter gives the Gauss–Newton Jacobian
of a low-rank parameterisation, and the adjoint and tangent agree to
rounding.

## Parameters

The gradient comes per row of the state and is summed into physical
parameters:

| Element | Parameters |
| --- | --- |
| Maxwell mantle cell | $\log\eta_M$ |
| Burgers mantle cell | $\log\eta_M$, $\log\tau_K = \log(\eta_K/\mu_K)$, $\log(\mu_K/G)$ |
| creeping fault patch | $\log\eta'$, and for a Burgers fault $\log\tau_K$, $\log(\mu_K/k_0)$ |

The Kelvin parameters are a relaxation time and a relative strength, not
$\eta_K$ and $\mu_K$, because postseismic data barely separate the latter.
The parameterisation of the field itself, for example a few
Karhunen–Loève modes of a smooth covariance, is left to the
caller, who applies it by the chain rule.

:::{admonition} In the code
:class: in-code

{mat:func}`solveLinearSequence` is the forward,
{mat:func}`computeAdjointGradient` the gradient per row and
{mat:func}`computeTangentSequence` the tangents;
{mat:func}`extractViscosityGradient` and {mat:func}`assembleRowPerturbation`
map rows to parameters and back. {mat:func}`assembleStationKernels` builds
$H$ and {mat:func}`computeDisplacementMisfit` the misfit, its residuals and
the Jacobian. They are in {doc}`../api/inverse`.
:::
