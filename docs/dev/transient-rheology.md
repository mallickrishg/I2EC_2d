# Transient (Burgers and Andrade) mantle rheology

Feature notes for branch `transient-rheology` (priority 1 of the feature
order: transient rheology → gravity → adjoint → nonlinear spin-up). Gravity and
the adjoint build on the interface below, so it is written to be extended.

## Status

- [x] Plan (this page)
- [x] Burgers in the semi-analytical solver (`assembleLinearSystem`)
- [x] Burgers in the numerical solver (`odeViscoelastic` and friends)
- [x] Andrade as a Prony series (`computeAndradePronySeries`)
- [x] Tests (`tests/unit/TransientRheologyTest.m`)
- [x] Example 5 and its docs page

## Rheology

Each mantle cell has a **Maxwell element** (the existing steady-state
dashpot, power law or Newtonian) in series with **p Kelvin elements** (a
spring μ_K in parallel with a dashpot η_K). With the code's convention for the
Maxwell dashpot (ε̇ = α σ, α = 1/η for n = 1; no factor 2), each Kelvin element
k obeys

    σ = μ_K,k ε_K,k + η_K,k ε̇_K,k        (both components 22 and 23)

and the strain rate that loads the kernels is the total

    ε̇ = ε̇_M + Σ_k ε̇_K,k,     ε̇_M = α τ^(n-1) σ.

| Rheology | Kelvin elements per cell | Solvers |
| --- | --- | --- |
| Maxwell (default) | 0 | both; results unchanged |
| Burgers | 1 | semi-analytical (n = 1), numerical (any n) |
| Andrade | p (Prony series, `computeAndradePronySeries`) | semi-analytical (n = 1), numerical (any n) |
| "power-law Burgers" | 1 or p, with n > 1 | numerical only |

Units are those of the rest of the code: η_K in MPa s (the same convention
as 1/α), μ_K in MPa (per unit strain, as the `LL` kernel).

### Selecting a rheology per cell

Two new properties of `geometry.shearZoneReceiver`, both [M × p]:

- `shz.etaK`: Kelvin viscosities (MPa s). `Inf` means "no element" in that
  cell, so p columns can describe cells with fewer elements.
- `shz.muK`: Kelvin moduli (MPa), positive where `etaK` is finite.

The constructor sets both to `zeros(M, 0)`: Maxwell everywhere. A cell is
Maxwell when every entry of its row is `Inf`; Burgers when one entry is
finite; Andrade when its row holds a Prony series. Mixed meshes (e.g. Andrade
wedge, Maxwell oceanic mantle) are rows with different content.
`extractKelvinElements(shz)` validates and returns the two arrays (also for
structs and for objects saved before the properties existed).

## Semi-analytical solver: augmented linear state

`sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change)` builds
dy/dt = A y + c and the coseismic jump b, which go unchanged into
`computePeriodicInitialCondition` and `solveAnalyticIvp`. The state is all
**rates**:

    y = [v (creeping patches); ε̇_M,22; ε̇_M,23;
         ε̇_K,22 and ε̇_K,23 of element 1 on the cells that have it; ...]

Every state entry maps to one "physical" entry of the Maxwell state
x = [v; ε̇_22; ε̇_23] (n0 = Nc + 2M), `sys.map`. The stressing rate depends
only on the total rates P y (P sums the entries that map to the same
physical entry, `sys.P`):

    dσ/dt = K (P y - y_pl) + σ̇_lock,   K = assembleStressKernel(evl, locked)

so, with η_row the viscosity of the element a row belongs to,

    dy_r/dt = (dσ/dt)_map(r) / η_row(r) - (μ_K/η_K)(r) y_r
    A = K(map, map) ./ η_row - diag(μ_K/η_K)
    c = σ̇_lock(map) ./ η_row - (K(map,:) ./ η_row) y_pl
    b = Δσ(map) ./ η_row

With no Kelvin elements map = 1:n0 and these are the expressions of
`examples/newtonianCycle.m` operation for operation, so Maxwell results are
bitwise unchanged.

- The Kelvin strain is recovered as ε_K = (σ - η_K ε̇_K)/μ_K with
  σ = η_M ε̇_M. It is continuous through the earthquake (Δσ/η_M·η_M -
  Δσ/η_K·η_K = 0), so in the periodic cycle the integral of ε̇_K over a cycle
  is zero: **the cycle-mean total strain rate and slip rate are those of the
  Maxwell model**, i.e. the long-term rates when the forcing balances.
- Stability: Maxwell and Kelvin elements are springs and dashpots, a passive
  (RC-like) network, so the spectrum is real and negative whenever the Maxwell
  operator is. μ_K > 0 is required: μ_K = 0 would be a dashpot in series, a
  neutral mode.
- Cells with infinite η_K carry no Kelvin state, so η_K → ∞ reduces to
  Maxwell exactly (and large finite η_K continuously).
- Cost: n = Nc + 2M(1 + p̄), with p̄ the mean number of elements per cell.
  Dense `eig` is O(n³), and computePeriodicInitialCondition and
  solveAnalyticIvp each do one: on the example mesh Maxwell n = 1061 (1 s),
  Burgers n = 2103 (3 s), 7-term Andrade n = 8355 (110 s). Keep Prony series short in
  the semi-analytical solver, or use them only where they matter.

For gravity: a restoring term adds to K (or to A via its own rows); the map
construction is unaffected. For the adjoint: A, P, map and η_row are all that
is needed; A is real with a real spectrum, Aᵀ has the same modes.

## Numerical solver: Kelvin spring-stress state

The ODE state of `odeViscoelastic` gains a block after the mantle stresses,
the spring stress s_K = μ_K ε_K of each element (MPa, the scale of the mantle
stresses, so the one ode45 tolerance suits both; storing ε_K itself, ~1e-6,
under AbsTol = 1e-10 lost accuracy):

    Y = [fault state (N·dgf); s22, s23 interleaved (M·dgf);
         sK22(:,1); sK23(:,1); sK22(:,2); sK23(:,2); ...]   (2·M·p)
    dε_K/dt = (s - s_K)/η_K,   ds_K/dt = μ_K dε_K/dt

(every cell, for each of the p columns; cells with η_K = Inf have zero rate).
Indexing of the existing entries is unchanged, so stress steps
(`integrateImposedCycles`) and all Maxwell code paths are untouched; with
p = 0 the block is empty and the arithmetic is identical. The Maxwell element
keeps the power law, so "power-law Burgers/Andrade" works here.

- `runImposedEarthquakeCycles` starts the Kelvin elements relaxed,
  s_K = σ_0 (zero Kelvin rate in the long-term state), and returns the states
  as a fifth output.
- `extractRatesFromOdeState` returns the total strain rate (Maxwell + Kelvin)
  and, as a fourth output, the spring stresses.
- Short Kelvin times make the system stiff; `integrateImposedCycles` and
  `runImposedEarthquakeCycles` take an optional ODE solver (default
  `@ode45`, e.g. `@ode15s`).
- Stress errors at the ode45 tolerance (1e-10 MPa) give strain-rate errors of
  1e-10/η_K: with η_K = 1e17 Pa s that is 1e-21 1/s, so the long-term state
  is kept to ~1e-6 of the rates instead of ~1e-11 for Maxwell alone.

## Andrade as a Prony series

Andrade transient creep under constant stress σ is ε_A(t) = β σ t^(1/3)
(β in MPa⁻¹ s^(-1/3)). A sum of Kelvin elements creeps as
Σ J_k (1 - exp(-t/τ_k)) σ, J_k = 1/μ_K,k, τ_k = η_K,k/μ_K,k.
`computeAndradePronySeries(beta, tmin, tmax)` takes τ_k log-spaced from
half a decade below tmin to half a decade above tmax ('Margin', default 0.5;
'TermsPerDecade', default 1) and fits J_k ≥ 0 by non-negative least squares to
the relative error of the creep function at 40 points per decade on
[tmin, tmax]. Since t^(1/3) has a fixed shape, the fit is done once and scaled
by β, so per-cell β costs nothing.

Without the margin (τ_k only inside [tmin, tmax]) the ends are fitted badly:
10% error at tmin and 1.5% in the interior whatever the number of terms. The
margin fixes both. Max relative error of the creep function on [tmin, tmax]
(`fit.maxRelErr`):

| Range | Terms per decade | Margin (decades) | Terms | Max rel. error |
| --- | --- | --- | --- | --- |
| 1e-3 – 2500 yr | 1 | 0.5 | 9 | 5.8e-3 |
| 0.1 – 1000 yr | 1 | 0.5 | 6 | 8.3e-3 |
| 1e-3 – 2500 yr | 2 | 0.5 | 16 | 4.0e-3 |
| 0.1 – 1000 yr | 2 | 0.5 | 11 | 4.1e-3 |
| 1e-3 – 2500 yr | 2 | 1 | 17 | 5.9e-5 |
| 0.1 – 1000 yr | 2 | 1 | 12 | 6.5e-5 |
| 1e-3 – 2500 yr | 3 | 1 | 24 | 2.3e-5 |
| 0.1 – 1000 yr | 3 | 1 | 17 | 1.0e-6 |

Example 5 uses 0.1 – 1500 yr with the defaults: 7 terms, 4.5e-3. Terms with a
zero fitted compliance are dropped. The tests check < 1e-2 (defaults) and
< 1e-4 (2 per decade, margin 1) on 1e-3 – 2500 yr.

Outside the range the series departs from Andrade by construction: below
tmin it creeps more slowly (t^(1/3) has an infinite initial rate), and above
tmax it saturates at Σ J_k (Andrade creep grows without bound). Choose tmax
to be a few cycle lengths so the whole cycle sees Andrade behaviour; the
saturation is what makes a periodic state exist.

## Results on the example mesh (example 5)

Cycle-mean departure from the long-term rates (slip, strain): Maxwell 4e-10,
1e-8; Burgers 2.5e-9, 2e-8; Andrade 2e-10, 8e-10. Largest Re(λ): -8.5e-12
(Maxwell, Burgers), -6.1e-12 1/s (Andrade). Mean mantle strain rate 0.01 yr
after the earthquake, over long-term: 11 (Maxwell), 114 (Burgers), 978
(Andrade).

## Verification (`tests/unit/TransientRheologyTest.m`)

- Maxwell default: `assembleLinearSystem` reproduces the operator, forcing and
  jump of the Maxwell solver bitwise (also checked on the example mesh:
  `newtonianCycle` before and after the change gives identical A, c, b, y0
  and solutions); the numerical RHS is unchanged.
- η_K = Inf gives the Maxwell system exactly; η_K → large converges to it.
- Single cell, constant stress (creep) and fixed total strain (relaxation):
  both solvers match the closed-form Burgers responses.
- Andrade single-cell creep matches β σ t^(1/3) within the fit error.
- Burgers and Andrade on the small model: no growing modes, real spectrum;
  cycle-mean strain and slip rates equal the long-term rates (with the slip
  deficit of a cycle released); the semi-analytical periodic cycle matches the
  numerical spun-up cycle (ode45, slip and strain rates to 1e-6).
- Power-law (n = 3) Burgers with two Kelvin columns starts and stays in the
  long-term state (ode45 and ode15s).
- `ExampleScriptsTest/transientRheologyExampleRuns` (slow tier) runs example 5
  and checks its cycle-mean departures are < 1e-6.

## Open / later

- Kelvin elements with a power law (nonlinear transient creep) — numerical
  only; not needed yet.
- A cheaper semi-analytical path for long Prony series (the Kelvin blocks
  are diagonal: a Schur-complement or rational eigen-solver could avoid the
  full O(n³) `eig`).
