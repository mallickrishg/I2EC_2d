# Example 5: transient mantle rheology (Burgers, Andrade)

The cycle of {doc}`example2` with three mantle rheologies that share the same
steady-state viscosity, $10^{19}$ Pa s, and so the same long-term flow:

| Rheology | Mantle cell | Parameters |
| --- | --- | --- |
| Maxwell | steady-state dashpot $\eta_M$ | as {doc}`example2` |
| Burgers | plus one Kelvin element | $\eta_K = 10^{18}$ Pa s, $\mu_K = 30$ GPa |
| Andrade | plus transient creep $\beta\sigma t^{1/3}$ (Prony series) | $\beta$: transient strain = elastic strain after 1 yr |

The equations stay linear, so the periodic cycle is still solved directly,
without spin-up (see {doc}`../theory/rheology`). Run it with

```matlab
run examples/example5_transientRheology
```

It takes about 10 minutes, nearly all of it for Andrade (see step 3).

## 1. Kelvin elements per cell

Transient creep is set on the mantle mesh by two $M \times p$ arrays, one
column per Kelvin element: viscosities `shz.etaK` (MPa s, `Inf` for no
element) and moduli `shz.muK` (MPa). Both are empty by default, which is the
Maxwell model. Cells can differ, e.g. Andrade in the wedge and Maxwell
elsewhere.

```matlab
% Burgers: one Kelvin element per cell
burgers.etaK = 1e18*1e-6 * ones(shz.N, 1);
burgers.muK  = m.mu * ones(shz.N, 1);

% Andrade: a Prony series fitted to beta*t^(1/3) on [0.1 yr, 3 cycles]
beta = 1/(m.mu*(1*yr)^(1/3));
fit = computeAndradePronySeries(beta*ones(shz.N, 1), 0.1*yr, 3*Trecur);
andrade.etaK = fit.etaK;
andrade.muK  = fit.muK;
```

{mat:func}`computeAndradePronySeries` places relaxation times one per decade
from half a decade below to half a decade above the fit range, and fits their
compliances by non-negative least squares:

```text
Andrade: 7 Kelvin elements, tau = 0.032 to 4.7e+03 yr, fit error 4.5e-03 on [0.1, 1500] yr
```

```{figure} ../_static/examples/example5_transientRheology_1.png
:alt: Creep curves of one cell for Maxwell, Burgers and Andrade, with the Prony fit
:width: 80%

Creep of one cell under constant stress, as viscous strain over elastic
strain. Burgers adds a transient that saturates after a few years, Andrade one
that grows as $t^{1/3}$. The Prony series (dashed) follows Andrade to 0.5%
within the fit range (dotted lines); below it the series creeps too slowly,
and beyond it the series saturates.
```

## 2. Solve the three cycles

`newtonianCycle` keeps the Kelvin elements it finds on `shz`, and builds the
linear system with {mat:func}`assembleLinearSystem`:

```matlab
for r = 1:3
    shz.etaK = rheo{r,2}.etaK; shz.muK = rheo{r,2}.muK;
    cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);
    sol = cyc.solve(tvec);           % total strain rates: Maxwell + Kelvin
end
```

`assembleLinearSystem` appends the Kelvin strain rates to the state: each
row is divided by the viscosity of its own element, the kernels act on the
total strain rate, and Kelvin rows also relax at $\mu_K/\eta_K$. The operator,
forcing and coseismic jump then go unchanged into
{mat:func}`computePeriodicInitialCondition` and {mat:func}`solveAnalyticIvp`.
`sol` holds the total strain rates; `extractRatesFromLinearState(Y, cyc.sys)`
also returns each element's share.

## 3. The cycle mean is the long-term rate

The Kelvin elements change *when* the mantle deforms, not how much: their
strain is continuous through the earthquake, so over a periodic cycle their
strain rates integrate to zero. The script checks this from the integrated
slip and strain over a whole cycle:

```text
Maxwell   1061 states, cycle-mean departure from the long-term rates: slip 2.7e-13, strain 1.3e-12  (2 s)
Burgers   2103 states, cycle-mean departure from the long-term rates: slip 2.3e-13, strain 3.4e-12  (10 s)
Andrade   8355 states, cycle-mean departure from the long-term rates: slip 2.9e-12, strain 9.9e-12  (491 s)
```

Each Kelvin element adds $2M$ states, and the eigendecomposition costs
$O(n^3)$: Burgers doubles the state, the 7-term Andrade series multiplies it
by eight. Keep Prony series short, or use them only in the cells that need
them.

## 4. Postseismic response

```{figure} ../_static/examples/example5_transientRheology_2.png
:alt: Mean mantle strain rate and mantle-driven surface velocity after the earthquake
:width: 100%

Left: area-weighted mean of the mantle strain-rate magnitude, over its
long-term value. Right: horizontal surface velocity at x = 200 km from mantle
flow alone (the fault's afterslip, faster and larger, is left out).
```

Just after the earthquake the transient elements make the mantle respond 10×
(Burgers) to 100× (Andrade) faster than the Maxwell mantle:

```text
mean mantle strain rate / long-term at 0.01 yr: [10.5 113 967]; at 500 yr: [1.3 1.3 1.29]
```

Burgers releases its transient within a few years (its Kelvin time is about
a year) and then creeps more slowly than the Maxwell mantle, which still has
its own relaxation to do, until the curves meet after about a century.
Andrade has no single time scale and relaxes steadily over decades. Late in the cycle the three agree, as they must, since they share the
steady-state viscosity and long-term loading.

## Numerical solver

The same Kelvin elements work in the numerical solver
({mat:func}`runImposedEarthquakeCycles`, {mat:func}`odeViscoelastic`), where
the Maxwell element may also be a power law ($n > 1$). The ODE state gains the
spring stresses $\mu_K\epsilon^K$ after the mantle stresses, the Kelvin
elements start relaxed, and short relaxation times make the system stiff, so
both drivers take an optional solver:

```matlab
shz.n = 3*ones(shz.N, 1);              % power-law Burgers
shz.etaK = burgers.etaK; shz.muK = burgers.muK;
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, ...
    stress_change, nreps, Trecur, @ode15s);
```

For a linear mantle the two solvers agree: the tests
`TransientRheologyTest/*PeriodicCycleMatchesSpunUpNumericalCycle` compare the
semi-analytical periodic cycle with the numerical spun-up cycle for Burgers
and Andrade.
