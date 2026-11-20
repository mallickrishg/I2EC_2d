# Example 1: power-law mantle, numerical solver

```matlab
run examples/example1_numericalPowerLaw
```

The mantle is a single power-law material,
$\dot\varepsilon = \alpha\,\tau^{n-1} s$ with $n = 3$ and uniform $\alpha$,
and the creeping part of the megathrust is rate-strengthening,
$d\log v/dt = \dot\tau/(a-b)\sigma$ with $(a-b)\sigma = 0.2$ MPa. The
geometry, locking and earthquake are those of the shared model (see
{doc}`index`). Because $\alpha$ is uniform, the $n = 3$ corner flow is the
long-term solution of this mantle; $\alpha$ is chosen so that the effective
viscosity is $10^{19}$ Pa s at the median long-term strain rate. The long-term
effective viscosity then ranges from $3\times10^{17}$ to $7\times10^{19}$ Pa s.
The coseismic stress step on the creeping fault is capped at 5 MPa next to
the locked zone.

The equations are nonlinear, so {mat:func}`runImposedEarthquakeCycles`
integrates them in time with `ode45`, from the long-term state through
15 cycles, and the last cycle is plotted. The script prints the change over
the last cycle as a check of periodicity: at most $3.4\times10^{-5}$ per cell
(median $2\times10^{-7}$).

```{figure} ../_static/examples/example1_numericalPowerLaw_1.png
:alt: Slip rate on the megathrust over the last cycle
:width: 80%

Slip rate on the megathrust over the last cycle (log time axis). Below the
locked zone ($x > 120$ km) afterslip starts at up to 100 times the plate rate
and decays over the cycle.
```

```{figure} ../_static/examples/example1_numericalPowerLaw_2.png
:alt: Mantle strain rate long-term and at three times in the last cycle
:width: 100%

Mantle strain rate: long-term (top left) and at three times in the last
cycle. With $n = 3$ the stressed regions weaken, so the postseismic strain
rates decay faster than in the Newtonian {doc}`example2`: by $t/T = 0.1$
they are close to their late-interseismic level.
```

```{figure} ../_static/examples/example1_numericalPowerLaw_3.png
:alt: Surface velocity at three times in the last cycle
:width: 80%

Horizontal (top) and vertical (bottom) surface velocity in the last cycle,
normalised by the plate rate.
```
