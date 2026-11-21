# Example 2: Newtonian mantle, semi-analytical solver

The model of {doc}`example1` with a Newtonian mantle ($10^{19}$ Pa s) and a
linear viscous creeping fault. The equations are then linear,
$dy/dt = Ay + c$, so the periodic cycle is solved directly, with no spin-up,
and evaluated exactly at any time. Run it with

```matlab
run examples/example2_analyticNewtonian
```

The steps below are the script and its helper `examples/newtonianCycle.m`,
shortened.

## 1. Load the model

```matlab
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; Vpl = m.Vpl; yr = m.yr;
```

As in {doc}`example1`: the megathrust, the mantle mesh and the cached kernels.

## 2. Set the rheology and solve the periodic cycle

```matlab
cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);
```

`newtonianCycle` sets the rheology on `rcv` and `shz` and builds the linear
system. Its core:

```matlab
% rheology: linear viscous fault, Newtonian mantle, locked above 30 km
rcv.faultLaw = 'linear';
rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) .* ones(rcv.N, 1);  % eta' (MPa s/m)
shz.n = ones(shz.N, 1);
shz.alpha = 1/(1e19*1e-6) .* ones(shz.N, 1);            % 1/viscosity (1/(MPa s))

% dy/dt = A*y + c, with y = [v on creeping patches; e22dot; e23dot]
eta = [rcv.Asigma(~locked); 1./shz.alpha; 1./shz.alpha];
A = assembleStressKernel(evl, locked) ./ eta;
c = c_locking - A*[rcv.Vpl(~locked); shz.e22pl; shz.e23pl];
b = [dtau(~locked); dsigma22; dsigma23] ./ eta;          % jump at the earthquake

% periodic cycle: the state after the earthquake repeats every Trecur
y0 = computePeriodicInitialCondition(A, c, b, Trecur);
```

{mat:func}`assembleStressKernel` stacks the fault and mantle kernels with the
locked patches removed; `c_locking` is the stressing from the locked zone and
the long-term loading enters through `c`. The script prints the decay times of
$A$'s slowest and fastest modes:

```text
slowest mode decays over 3715 yr, fastest over 0.008 yr (cycle: 500 yr)
```

## 3. Evaluate the cycle at any time

```matlab
tfrac = [0.01 0.1 0.99];
sol = cyc.solve(tfrac*cyc.Trecur);
```

`cyc.solve` calls {mat:func}`solveAnalyticIvp`, which evaluates the solution
exactly at the requested times (one matrix–vector product each), and splits it
with {mat:func}`extractRatesFromLinearState`. `sol.V` is the slip rate
(patches × times, 0 on the locked zone), `sol.e22dot` and `sol.e23dot` the
mantle strain rates (cells × times); `sol.slip`, `sol.e22`, `sol.e23` are
their integrals since the earthquake.

## 4. Plot the results

**Mantle strain rate**, long-term and at three times, drawn with
`plotShz2d(shz, edot(:,i))`. Just after the earthquake the whole mantle deforms
faster than in the long term; the excess decays through the cycle.

```{figure} ../_static/examples/example2_analyticNewtonian_1.png
:alt: Mantle strain rate long-term and at three times after the earthquake
:width: 100%
```

**Slip rate** on the creeping part of the megathrust, `sol.V(~locked,:)`.
Afterslip runs above the plate rate for the first decades.

```{figure} ../_static/examples/example2_analyticNewtonian_2.png
:alt: Slip rate on the creeping megathrust at three times
:width: 80%
```

**Surface velocity**, the same kernel product as in {doc}`example1`:

```matlab
G = @(j) [m.devl.KO(:,:,j), m.devl.LO(:,:,j,1), m.devl.LO(:,:,j,2)];
dstate = [sol.V; sol.e22dot; sol.e23dot] - [rcv.Vpl; shz.e22pl; shz.e23pl];
vx = G(1)*dstate - m.Ghx*(m.hinge.Vpl*Vpl);     % j = 2 for vz
```

```{figure} ../_static/examples/example2_analyticNewtonian_3.png
:alt: Horizontal and vertical surface velocity at three times
:width: 80%
```
