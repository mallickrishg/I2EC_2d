# Rheology

The stress sets how fast the fault slips and the mantle flows. This page
gives the laws for the mantle cells and the fault patches, and how locking
and earthquakes enter.

## Mantle

Each cell flows as a power law in the deviatoric stress $\tau$:

$$
\dot{\epsilon}_{ij} = \alpha\, \tau^{\,n-1}\, \sigma_{ij},
\qquad \tau = \sqrt{\sigma_{22}^2 + \sigma_{23}^2}.
$$

With $n = 1$ the mantle is Newtonian and $\alpha = 1/\eta$.

## Transient creep (Burgers, Andrade)

Each cell may also have Kelvin
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

## Fault

The creeping fault follows one of two laws:

| Law | Equation | Used by |
| --- | --- | --- |
| rate-strengthening friction | $\dfrac{d \log v}{dt} = \dfrac{\dot{\tau}}{(a-b)\sigma}$ | numerical solver (default) |
| linear viscous | $\dfrac{dv}{dt} = \dfrac{\dot\tau}{\eta'}$ | both solvers |

## Burgers fault (transient afterslip)

A creeping patch may also have Kelvin
elements in series with its dashpot, exactly as a mantle cell: the slip of
each obeys $\tau = \mu_{K,k} s^K_k + \eta_{K,k} \dot s^K_k$, and the slip rate
that loads the kernels is $v + \sum_k \dot s^K_k$. The Kelvin slip is
continuous through an earthquake, so over a periodic cycle the Kelvin slip
rates average to zero and the long-term slip rate is unchanged.

## Locking and earthquakes

Locked patches do not slip between earthquakes.
Every recurrence interval $T$ an earthquake releases their slip deficit,
$v^{\infty} T$, and its stress change is applied as an instantaneous step.

:::{admonition} In the code
:class: in-code

the mantle takes `shz.alpha` and `shz.n`, and Kelvin elements in
`shz.etaK` (MPa s; `Inf` for none) and `shz.muK` (MPa), $M \times p$ each
(default $M \times 0$: Maxwell). {mat:func}`computeAndradePronySeries`
returns them for Andrade creep; {doc}`../rheology` shows how to set each
rheology and call the solvers, and {doc}`../examples/example5` compares the three.
The fault takes
`rcv.faultLaw` (`'ratestrengthening'` or `'linear'`, read by
{mat:func}`isLinearFaultLaw`) and `rcv.Asigma`, which holds $(a-b)\sigma$ (MPa)
or $\eta'$, a viscosity per unit length (MPa s/m). Fault Kelvin elements are
`rcv.etaK` (MPa s/m; `Inf` for none) and `rcv.muK` (MPa/m), $N \times p$
each (default $N \times 0$). Locked patches are `rcv.pinnedPosition`.
:::
