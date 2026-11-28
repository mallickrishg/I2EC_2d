# Solvers

The governing equations are a system of ordinary differential equations.
I2EC integrates them in time, or, when the model is linear, solves the
periodic cycle exactly.

| | Numerical | Semi-analytical |
| --- | --- | --- |
| Rheology | any (power law, friction; Burgers, Andrade) | linear: $n = 1$ (Maxwell, Burgers, Andrade), linear viscous fault |
| Method | time integration with `ode45` | eigendecomposition, exact in time |
| Periodic cycle | spin-up over many cycles | solved directly |
| Example | {doc}`../examples/example1` | {doc}`../examples/example2` |

{doc}`../examples/example4` runs both on the same linear model; they agree to
about $10^{-9}$.

## Numerical integration

The equations are integrated with MATLAB's adaptive Runge–Kutta solver
`ode45` (tolerance $10^{-10}$), restarting after each earthquake. Repeating the
cycle from the long-term state spins the model up until it is periodic.

:::{admonition} In the code
:class: in-code

{mat:func}`runImposedEarthquakeCycles` repeats the cycle;
{mat:func}`integrateImposedCycles` integrates one sequence of events and
{mat:func}`odeViscoelastic` is the right-hand side. Spin-up convergence is
not checked: compare the start and end of the last cycle, as example 1 does.
Kelvin elements add their spring stresses $\mu_K\epsilon^K$ (on the fault
$\mu_K s^K$) to the state;
short relaxation times make the system stiff, so both drivers take an
optional solver such as `@ode15s`.
:::

## Semi-analytical solution

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
kernels guarantee this ({ref}`mantle-kernels`).

:::{admonition} In the code
:class: in-code

{mat:func}`assembleStressKernel` builds $K$ without the
locked patches, and {mat:func}`assembleLinearSystem` builds $A$, $c$ and $b$.
With Kelvin elements it appends their strain rates (and fault Kelvin slip
rates) to $y$: each row is
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
