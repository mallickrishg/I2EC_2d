# Set up the rheology and solve

The rheology lives on two objects: the megathrust `rcv`
({mat:class}`geometry.receiver <receiver>`) has one value per fault patch, and
the mantle mesh `shz` ({mat:class}`geometry.shearZoneReceiver <shearZoneReceiver>`)
one value per triangular cell. Set their properties, then call a solver. This
page goes through the steps in order; the code blocks run as one script.

| Object | Property | Size | Units | Meaning |
| --- | --- | --- | --- | --- |
| `rcv` | `faultLaw` | | | `'ratestrengthening'` (default) or `'linear'` |
| `rcv` | `Asigma` | N × 1 | MPa, or MPa s/m | $(a-b)\sigma$, or the viscosity per unit length of a linear fault |
| `rcv` | `pinnedPosition` | N × 1 | logical | locked patches |
| `rcv` | `Vpl` | N × 1 | m/s | long-term slip rate |
| `shz` | `n` | M × 1 | | power-law exponent (1: linear) |
| `shz` | `alpha` | M × 1 | 1/(MPa$^n$ s) | fluidity; $1/\eta$ for $n = 1$ |
| `shz` | `etaK`, `muK` | M × p | MPa s, MPa | Kelvin elements (transient creep); $M \times 0$ for none |
| `shz` | `e22pl`, `e23pl` | M × 1 | 1/s | long-term strain rates |

Viscosities are in **MPa s**: multiply a value in Pa s by $10^{-6}$.

## 1. Geometry and kernels

```matlab
m = exampleModel('StressKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl; Vpl = m.Vpl; yr = m.yr;
Trecur = 500*yr;
```

`exampleModel` loads the example geometry and its cached kernels. For your own
geometry, read the files described in {doc}`inputs` and compute the kernels
once with {mat:func}`computeAllStressKernelsBem`.

## 2. The fault

```matlab
locked = abs(rcv.xc(:,2)) < 30e3;          % locked above 30 km depth
rcv.pinnedPosition = locked;
rcv.Vpl = Vpl .* ones(rcv.N, 1);

% linear viscous fault: eta' = 5e17 Pa s spread over the fault width
rcv.faultLaw = 'linear';
rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) .* ones(rcv.N, 1);
```

The semi-analytical solver needs the linear fault. The numerical solver also
takes rate-strengthening friction: `rcv.faultLaw = 'ratestrengthening'` with
`rcv.Asigma` $= (a-b)\sigma$ in MPa (see {doc}`examples/example1`).

## 3. The mantle

Every property is per cell, so it can vary through the mesh: use the cell
centres `shz.xc` (m, $z$ negative) to set it by region or depth.

**Maxwell (linear).** A viscosity per cell:

```matlab
eta = 1e19 * ones(shz.N, 1);                % Pa s
eta(shz.xc(:,2) < -100e3) = 1e20;          % stronger below 100 km
shz.n = ones(shz.N, 1);
shz.alpha = 1 ./ (1e-6*eta);               % 1/(MPa s)
shz.etaK = zeros(shz.N, 0);                % no transient creep
shz.muK  = zeros(shz.N, 0);
```

**Burgers.** One Kelvin element (spring $\mu_K$ and dashpot $\eta_K$ in
parallel) in series with the Maxwell dashpot, one column per element:

```matlab
shz.etaK = 1e-6*1e18 * ones(shz.N, 1);     % MPa s
shz.muK  = m.mu * ones(shz.N, 1);          % MPa
```

**Andrade.** Creep $\epsilon = \beta\sigma t^{1/3}$, approximated by several
Kelvin elements fitted over a range of times
({mat:func}`computeAndradePronySeries`, about 1% error):

```matlab
beta = 1/(m.mu*(1*yr)^(1/3));              % 1/(MPa s^(1/3))
fit = computeAndradePronySeries(beta*ones(shz.N, 1), 0.1*yr, 3*Trecur);
shz.etaK = fit.etaK;
shz.muK  = fit.muK;
```

Each Kelvin element adds two unknowns per cell, which sets the cost. On the
example mesh (521 cells) Maxwell has 1061 unknowns, Burgers 2103 and the
Andrade fit above (7 elements) 8355: its semi-analytical solve takes about
8 minutes and three numerical cycles about 10. Fewer terms
(`'TermsPerDecade'`) or a narrower time range make it cheaper. Set `etaK` to
`Inf` in cells that have no element.

**Power law.** Set `shz.n` and a fluidity per cell. Only the numerical solver
takes $n \ne 1$; {doc}`examples/example1` sets `alpha` from a reference
viscosity at the median long-term strain rate.

## 4. Long-term loading

```matlab
[e22, e23] = computeLongTermStrainRates(shz, rcv.dip(1)*pi/180, ...
                                        [0, 20e3], [-140e3, 35e3], 1);
shz.e22pl =  e22 .* Vpl;
shz.e23pl = -e23 .* Vpl;
```

{mat:func}`computeLongTermStrainRates` gives the corner-flow strain rates per
unit plate rate; its last argument is the power-law exponent (1 or 3). The two
translations shift the oceanic and wedge regions into the corner-flow frame
(see {doc}`theory`). The Kelvin elements do not change the long-term flow.

## 5. The earthquake

```matlab
slip = Trecur*Vpl .* locked;               % slip deficit of one cycle
stress_change = struct('Nevents', 1, 'Timing', 0, ...
    'dtau', evl.KK*slip, ...
    'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
```

## 6. Solve

| Rheology | Semi-analytical | Numerical |
| --- | --- | --- |
| Maxwell, Burgers, Andrade with $n = 1$ and a linear fault | yes, exact and fastest | yes |
| power law ($n \ne 1$) or rate-strengthening fault | no | yes |

### Semi-analytical solver

```matlab
sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change);
y0  = computePeriodicInitialCondition(sys.A, sys.c, sys.b, Trecur);
tvec = [0.01 0.1 0.99]*Trecur;             % times after the earthquake
[Y, U] = solveAnalyticIvp(sys.A, sys.c, y0, tvec);
[V, e22dot, e23dot, parts] = extractRatesFromLinearState(Y, sys);
```

{mat:func}`assembleLinearSystem` builds the linear system $dy/dt = Ay + c$ and
the jump $b$ of each earthquake, with the Kelvin elements of `shz` included;
{mat:func}`computePeriodicInitialCondition` finds the periodic cycle with no
spin-up, and {mat:func}`solveAnalyticIvp` evaluates it exactly at any times.
`V` is the slip rate (patches × times, 0 on locked patches), `e22dot` and
`e23dot` the total mantle strain rates (cells × times), and `parts` splits
them into the Maxwell and Kelvin elements. `U` holds the integrals since the
earthquake, split the same way:

```matlab
[slip_t, e22_t, e23_t] = extractRatesFromLinearState(U, sys);
```

### Numerical solver

```matlab
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, ...
                                stress_change, 3, Trecur, @ode15s);
```

{mat:func}`runImposedEarthquakeCycles` starts from the long-term state,
integrates the given number of cycles and returns the last one. It takes any
rheology on this page, including the power law and friction. Use `@ode15s`
when there are Kelvin elements (their short relaxation times make the system
stiff); the default `ode45` suits Maxwell. Run enough cycles to reach a
periodic state and check it by comparing the start and end of the last cycle,
as {doc}`examples/example1` does; {ref}`valid-range` says how many.

The examples put these steps together: {doc}`examples/example1` (power law and
friction, numerical), {doc}`examples/example2` (Maxwell, semi-analytical) and
{doc}`examples/example5` (Maxwell, Burgers and Andrade).
