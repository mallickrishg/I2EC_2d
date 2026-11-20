# Example 2: Newtonian mantle, semi-analytical solver

```matlab
run examples/example2_analyticNewtonian
```

The mantle is Newtonian with a uniform viscosity of $10^{19}$ Pa s and the
creeping fault is linear viscous, with the geometry, locking and earthquake of
the shared model (see {doc}`index`; the rheology is set by `newtonianCycle`).
The cycle is then a linear ODE, $dy/dt = Ay + c$, so its periodic state is
found directly with {mat:func}`computePeriodicInitialConditionBalanced`, with
no spin-up, and {mat:func}`analyticSolveIVP` evaluates the solution exactly at
$t/T = 0.01$, $0.1$ and $0.99$ after the earthquake. The script also prints
the slowest and fastest relaxation times of $A$: about 3700 years and 3 days.

```{figure} ../_static/examples/example2_analyticNewtonian_1.png
:alt: Mantle strain rate long-term and at three times after the earthquake
:width: 100%

Mantle strain rate: long-term (top left) and at three times after the
earthquake. Shortly after it the whole mantle deforms faster than in the long
term; the excess decays through the cycle.
```

```{figure} ../_static/examples/example2_analyticNewtonian_2.png
:alt: Slip rate on the megathrust at three times
:width: 80%

Slip rate on the creeping part of the megathrust (the locked zone, $v = 0$,
is not shown on the log scale).
```

```{figure} ../_static/examples/example2_analyticNewtonian_3.png
:alt: Horizontal and vertical surface velocity at three times
:width: 80%

Horizontal (top) and vertical (bottom) surface velocity, normalised by the
plate rate.
```
