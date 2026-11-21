# Example 1: power-law mantle, numerical solver

A subduction zone with a power-law ($n = 3$) mantle and a rate-strengthening
creeping fault, loaded by plate motion and an earthquake every 500 years. The
equations are nonlinear, so they are integrated in time. Run it with

```matlab
run examples/example1_numericalPowerLaw
```

The steps below are the script, shortened.

## 1. Load the model

```matlab
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl; Vpl = m.Vpl; yr = m.yr;
```

`exampleModel` (in `examples/`) adds the I2EC folders to the path and returns
the megathrust `rcv` (49 patches), the mantle mesh `shz` (521 triangles), the
plate rate and the kernels: `evl` for stress, `m.devl` for surface
displacement, loaded from `examples/kernels/`. To build kernels for your own geometry,
see {doc}`../inputs` and {mat:func}`computeAllStressKernelsBem`.

## 2. Set the rheology

The fault and the mantle are handle objects; set their properties directly.

```matlab
% fault: rate-strengthening, locked above 30 km
rcv.faultLaw = 'ratestrengthening';
rcv.Asigma = 0.2 .* ones(rcv.N, 1);            % (a-b)sigma, MPa
locked = abs(rcv.xc(:,2)) < 30e3;
rcv.pinnedPosition = locked;
rcv.Vpl = Vpl .* ones(rcv.N, 1);

% long-term loading: n = 3 corner flow
[e22_lt, e23_lt] = computeLongTermStrainRates(shz, rcv.dip(1)*pi/180, ...
                                          [0, 20e3], [-140e3, 35e3], 3);
shz.e22pl = e22_lt .* Vpl;
shz.e23pl = -e23_lt .* Vpl;

% mantle: power law edot = alpha * tau^(n-1) * s
shz.n = 3 .* ones(shz.N, 1);
eta_ref  = 1e19*1e-6;                           % MPa s
edot_ref = median(sqrt(shz.e22pl.^2 + shz.e23pl.^2));
shz.alpha = edot_ref ./ (eta_ref*edot_ref).^3 .* ones(shz.N, 1);
```

{mat:func}`computeLongTermStrainRates` gives the strain rates per unit plate rate.
`shz.alpha` makes the effective viscosity $10^{19}$ Pa s at the median
long-term strain rate; faster-deforming cells are weaker. The script prints the
resulting range:

```text
long-term effective viscosity: 2.7e+17 to 7e+19 Pa s (1e19 at 3e-16 1/s)
```

## 3. Impose the earthquake

Each earthquake releases the slip deficit of one cycle on the locked zone. Its
stress change on the creeping fault and in the mantle comes from the kernels:

```matlab
Trecur = 500*yr;
slip = Trecur*Vpl .* locked;
dtau = min(evl.KK*slip, 5);  dtau(locked) = 0;   % cap at 5 MPa
stress_change = struct('Nevents', 1, 'Timing', 1, 'dtau', dtau, ...
    'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
```

## 4. Integrate the cycles

```matlab
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, ...
                                                    stress_change, 15, Trecur);
```

{mat:func}`runImposedEarthquakeCycles` starts from the long-term state, runs
15 cycles with `ode45` (about a minute), and returns the last one: times `t` (s), slip rate `V`
(times × patches, m/s) and mantle strain rates `e22dot`, `e23dot` (times ×
cells, 1/s). Fifteen cycles bring it close to periodic; the script checks by
comparing the start and end of the last cycle:

```text
change over the last cycle: strain rate median 2.2e-07, max 3.4e-05; slip rate max 3.7e-09 (relative, per cell/patch)
```

## 5. Plot the results

**Slip rate** at three creeping patches, `V(:,ip)/Vpl` against `t/Trecur`.
Afterslip starts at up to 100 times the plate rate and decays below it within
the first tenth of the cycle; the shallowest patch, next to the locked zone,
ends slowest.

```{figure} ../_static/examples/example1_numericalPowerLaw_1.png
:alt: Slip rate at three creeping patches over the last cycle
:width: 80%
```

**Mantle strain rate**, $\sqrt{\dot\varepsilon_{22}^2 + \dot\varepsilon_{23}^2}$
at three times, drawn with `plotShz2d(shz, edot(i,:)')` and
`plotPatches2d(rcv)`. The stressed regions weaken, so the postseismic strain
rates decay faster than in the Newtonian {doc}`example2`.

```{figure} ../_static/examples/example1_numericalPowerLaw_2.png
:alt: Mantle strain rate long-term and at three times in the last cycle
:width: 100%
```

**Surface velocity** is a kernel product with the departure from the long-term
state, plus the kinematic hinge fault:

```matlab
G = @(j) [m.devl.KO(:,:,j), m.devl.LO(:,:,j,1), m.devl.LO(:,:,j,2)];
dstate = [V(i,:)'; e22dot(i,:)'; e23dot(i,:)'] - [rcv.Vpl; shz.e22pl; shz.e23pl];
vx = G(1)*dstate - m.Ghx*(m.hinge.Vpl*Vpl);     % j = 2 for vz
```

Five years after the earthquake the surface velocity alternates along the
profile. This comes from the power-law mantle (the Newtonian {doc}`example2`
has none of it): the same pattern, with amplitudes within about 0.6 of the
plate rate, appears on the meshes of {doc}`example3` with 495 and 1955 cells.

```{figure} ../_static/examples/example1_numericalPowerLaw_3.png
:alt: Surface velocity at three times in the last cycle
:width: 80%
```
