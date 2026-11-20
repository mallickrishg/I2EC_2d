# Example 4: semi-analytical against numerical

```matlab
run examples/example4_analyticVsNumerical
```

The script solves the cycle of {doc}`example2` with both solvers. The
numerical run uses the linear viscous fault law (`rcv.faultLaw = 'linear'`,
set by `newtonianCycle`), so both solve the same linear ODE. It starts from
the analytic state just before the earthquake and integrates one cycle with
{mat:func}`viscoeq_sequence_simulation` and {mat:func}`ode_viscoelastic`. It
must follow the analytic solution at every time and end where it started. On
the example they agree to $3.8\times10^{-9}$ in slip rate and
$3.6\times10^{-11}$ in strain rate (relative), and the numerical end state
matches its start to $3\times10^{-11}$.

```{figure} ../_static/examples/example4_analyticVsNumerical_1.png
:alt: Time series of fault slip rate and mantle strain rate from both solvers
:width: 100%

Slip rate at three creeping patches (32, 40 and 48 km depth) and strain rate
at three mantle cells through one cycle: semi-analytical (coloured) and
numerical (black dashed).
```

Starting the numerical run from the long-term state instead would need a
spin-up. The slowest mode of this model decays over about 3700 years, so the
remaining error shrinks only by a factor of 0.87 per 500-year cycle, and about
100 cycles are needed to reach $10^{-6}$ (printed by the script).
`TimeIntegrationTest/spunUp*` checks the spin-up on a smaller model.
