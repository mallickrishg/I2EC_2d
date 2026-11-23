# Example 4: semi-analytical against numerical

The cycle of {doc}`example2` solved twice: semi-analytically, and by
integrating the same equations in time with `ode45`, the solver of
{doc}`example1`. With a linear viscous fault (`rcv.faultLaw = 'linear'`) both
solve the same linear ODE, so they must agree at every time. Run it with

```matlab
run examples/example4_analyticVsNumerical
```

The steps below are the script, shortened.

## 1. Solve the cycle semi-analytically

```matlab
m = exampleModel('StressKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl;
cyc = newtonianCycle(rcv, shz, evl, m.Vpl, m.yr);
```

As in {doc}`example2`. The script also prints how slowly a numerical spin-up
from the long-term state would converge, from the slowest mode of `cyc.A`:

```text
slowest mode: 3715 yr (rho = 0.874 per cycle, 103 cycles to spin up to 1e-6)
```

So the numerical run below starts from the periodic state instead.

## 2. Integrate one cycle numerically

The numerical solver's state is slip rate and mantle *stress*, so convert the
analytic state just before the earthquake:

```matlab
locked = cyc.locked; Nc = sum(~locked); M = shz.N;
ypre = cyc.y0 - cyc.b;                          % [v(creeping); e22dot; e23dot]
Y0 = zeros(rcv.N + 2*M, 1);                     % [v; s22, s23 interleaved]
Y0(~locked) = ypre(1:Nc);
Y0(rcv.N+1:2:end) = ypre(Nc+1:Nc+M) ./ shz.alpha;
Y0(rcv.N+2:2:end) = ypre(Nc+M+1:end) ./ shz.alpha;

[t, Y] = integrateImposedCycles(@(t,y) odeViscoelastic(t,y,rcv,shz,evl), ...
    Y0, rcv, shz, cyc.stress_change, cyc.Trecur);
[Vn, e22n, e23n] = extractRatesFromOdeState(Y, rcv, shz);
```

{mat:func}`integrateImposedCycles` integrates
{mat:func}`odeViscoelastic` with `ode45`, applying the earthquake in
`cyc.stress_change` at $t = 0$, and
{mat:func}`extractRatesFromOdeState` turns the state back into slip
and strain rates.

## 3. Compare

```matlab
sol = cyc.solve(t');                            % analytic, at the same times
```

The script prints the largest difference over the cycle and whether the
numerical run ends where it started:

```text
largest difference over the cycle: slip rate 3.6e-09, strain rate 2.9e-11 (relative)
numerical end state vs start state: 7.1e-11 (relative)
```

## 4. Plot the results

Slip rate at three creeping patches and strain rate at three mantle cells
through the cycle: semi-analytical (coloured) and numerical (black dashed).
The curves lie on top of each other.

```{figure} ../_static/examples/example4_analyticVsNumerical_1.png
:alt: Time series of fault slip rate and mantle strain rate from both solvers
:width: 100%
```
