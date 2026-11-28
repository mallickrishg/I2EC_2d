# Adjoint gradients for viscosity inversion (linear rheologies)

Feature notes for branch `adjoint-inverse` (feature 3 in `tasks.md`; moved
ahead of gravity on 2026-11-25). Scope: Maxwell and Burgers
mantles, n = 1. Andrade is skipped for now (user, 2026-11-25): no β
parameter, no Andrade tests or example; the row-level core is written for any
number of Kelvin columns, so it can be added later. No power-law adjoint. No gravity yet, but gravity-corrected
kernels K_g must go through unchanged (see "Gravity" below).

## Status

This page is the feature's task list and record: plan, checklist, results
and open items all live here, on the branch. `tasks.md` keeps only its
one-line entry in the Features list, updated once when the branch merges
(user, 2026-11-25).

- [x] Plan (this page), derivation checked numerically on a random 6-state
      system (adjoint vs central FD, rel. error 1e-9; scratch prototype only)
- [x] User approval of the API and open questions (2026-11-25)
- [x] Sequence forward solve (`solveLinearSequence`)
- [x] Core adjoint VJP (`computeAdjointGradient`) and tangent
      (`computeTangentSequence`)
- [x] Row ↔ physical-parameter maps (`extractViscosityGradient`,
      `assembleRowPerturbation`)
- [x] Observation layer (`assembleStationKernels`, `computeDisplacementMisfit`)
- [x] Tests (`tests/unit/AdjointTest.m`)
- [x] Example 6 (synthetic Gauss–Newton inversion in KL modes) and docs page
- [x] `inverse/` folder wired into the test fixture, `exampleModel`, Sphinx and the docs workflow

## Done since the first merge-ready commit (pushed to `adjoint-inverse`)

- [x] `89b884f`: balanced eigenvectors without inverting V (accurate when
      viscosity spans decades); tangent directions in chunks of 12.
- [x] `2b8cd15`: separable pairwise time integrals
      (`computePairIntegralSums`): adjoint and tangent cost no longer
      depends on the time sampling (80 log-spaced epochs: 0.86 s, was 61 s;
      log-spaced converges to dense uniform at second order).
- [x] `dabe09e`: missing data (NaN / zero weight) in
      `computeDisplacementMisfit`; tests for GNSS-like messy sampling and
      stations switching on/off; example-mesh stress test (119 epochs,
      69 distinct gaps, 27% missing) matches FD to 7e-7.

## Next (start 2026-11-26): rework example 6

Uncommitted in the working tree: `examples/example6_adjointInversion.m`
(thermal truth, SVD resolution), its docs page and figures, README/docs
tables, `ExampleScriptsTest` and the example-6 notes on this page. None of
it is pushed (user's request).

- [ ] Truth as now: thermal model (20 Myr plate, wedge, Arrhenius with
      V* = 12 cm^3/mol: high-low-high eta_M, 1e18.9-1e22 Pa s),
      eta_K = eta_M/10, one mu_K = G/2. Recurrence 500 yr.
- [ ] Parameterisation with different scales per field (user): most
      modes for log tauK, fewer for log eta_M, a single value (0th mode)
      for log(mu_K/G). Isotropic Matern KL modes (least-biased prior).
- [ ] Data: GNSS-like daily series, log-downsampled, with messy gaps and
      stations that start late / stop early (NaN in obs.d).
- [ ] Resolution analysis (prior-preconditioned whitened Jacobian, SVD,
      lambda > 1, effective parameters), split of resolved directions
      between the fields; leading directions on the mesh.
- [ ] Keep it cheap (the 37-min and 45-min runs were too long); show the
      stations and use axis equal on mesh plots.
- [ ] Known issues to keep in mind: with tauK as a free field the
      optimizer drove tauK to seconds (below the first epoch, invisible
      to data) in an earlier run; Levenberg-Marquardt now uses Marquardt
      scaling.
- [ ] Possibly add a reusable resolution tool to inverse/ (full prior
      square root by default, matrix-free randomized SVD via tangent
      J*v and adjoint J'*u for big meshes).
- [ ] Make the adjoint the workhorse of the example (user, 2026-11-25:
      the explicit-Jacobian Gauss-Newton used so far never calls
      computeAdjointGradient). Decided: fminunc with the adjoint gradient
      ('SpecifyObjectiveGradient', true; quasi-Newton / 'lbfgs'): one
      forward + one adjoint per iteration at any number of parameters,
      including mesh-scale. Objective = data misfit + prior term, gradient
      = Phi' * (adjoint gradient per cell) + prior gradient. Not
      lsqnonlin. Compare cost and convergence with the explicit-Jacobian
      Gauss-Newton only as a reference.
- [ ] State the regularisation explicitly (user, 2026-11-25): the unit
      Gaussian prior on the KL coefficients is Tikhonov with the
      Matern precision; its strength sits in the prior std (sigma =
      1 decade for eta_M and tauK, 0.5 for mu_K/G), the length 60 km and
      nu = 3/2, all chosen ad hoc so far. Show the sensitivity to sigma,
      and/or choose it from the data (discrepancy chi^2/N ~ 1, L-curve, or
      the marginal likelihood). Note the LM damping is step control only.
- [ ] Then: update docs page, slow test, this page; commit and push
      example 6 only when the user says so.

## On hold (user, 2026-11-26; pick up when the user is back)

Branch plan (user, 2026-11-26): build the Burgers fault (forward model,
then its adjoint) on `adjoint-inverse`, test it and merge to main. Then
make a new branch `example6` for the example (fault inversion, tweaks,
docs page, slow test). Nothing is built yet.

- [ ] Burgers fault: a Kelvin element in series with the creeping-fault
      dashpot (today the fault is one dashpot, eta' per patch, no Kelvin
      rows). New fault Kelvin rows in assembleLinearSystem and the
      numerical solver, parameterised as the mantle (log eta', log tauK,
      log muK/mu), with tests. The adjoint is row-wise, so it mostly
      carries over; extend extractViscosityGradient /
      assembleRowPerturbation to the fault Kelvin rows.
- [ ] Then invert uniform fault Burgers parameters on the creeping section
      in example 6 (afterslip vs mantle trade-off); coseismic slip stays
      known (uniform slip deficit on the locked zone, 0-30 km).
- [ ] Example 6 tweaks still open: fits as data minus known terms, drop
      the all-cell rms error, cap iterations (~12 min now).
- [ ] Rename examples/newtonianCycle.m (it also sets up Burgers/Andrade
      cycles; suggestion linearCycle). About 15 files; rename on main
      after the merge, not on this branch. Undecided.

## Results (2026-11-25)

- `tests/unit/AdjointTest` (small model, n = 55 Maxwell / 103 Burgers):
  every per-parameter adjoint gradient (fault η′, η_M, and for Burgers τ_K,
  μ_K/G or the raw logs) matches central differences to ≤ 1e-6 of the
  largest gradient component (h = 1e-4; for weakly sensitive parameters the
  difference is FD round-off: it grows as 1/h); dot-product test
  tangent/adjoint to 1e-10; tangent vs FD 1e-6; modal periodic state vs
  `computePeriodicInitialCondition` 1e-9; divided differences vs Opitz
  `expm` 1e-13, including equal, near-equal, zero and complex nodes.
- Example 6 (example mesh, Burgers, n = 2103, 21 stations × 2 components ×
  80 epochs): all three fields vary in space, parameterised as log η_M,
  log τ_K, log μ_K/G with 5, 4, 3 Matérn KL modes. Case A inverts all three
  (12 coefficients), case B holds μ_K fixed and inverts η_M and η_K through
  τ_K (9). Forward + 12 tangents ≈ 5 s per iteration; adjoint gradient
  equals Jᵀr from the tangents to 1e-15; χ²/N = 1.02 after 3 Gauss–Newton
  iterations in both cases, |ξ − ξ_true| = 0.08 (A) and 0.05 (B). (The
  first version, 2026-11-25, was a Maxwell-only inversion; the user asked
  for Burgers with spatially variable η_M, η_K, μ_K and the fixed-μ_K case.)
- Not done: Andrade (deferred by the user); the adjoint is written for any
  number of Kelvin columns but tested with one.

## 1. Forward model

`sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change)` gives,
for every rheology,

    dy/dt = A y + c,   A = D K̃,  c = D f,  jump at an event: y⁺ = y⁻ + b,  b = D Δσ
    D = diag(1/η_row),  K̃ = K(map,map) - diag(μ_K on Kelvin rows, 0 elsewhere)

(`sys.eta` is η_row; `sys.decay` = μ_K/η_K). Every row of A, c and b is
divided by the viscosity of the element it belongs to — fault patch
(`rcv.Asigma`), Maxwell dashpot (1/`shz.alpha`, two rows per cell) or
Kelvin dashpot (`shz.etaK`, two rows per cell and element). Nothing else
depends on η: f (locked-zone loading minus K y_pl), Δσ (event stress
changes) and K̃ are fixed.

**Time line.** Two parts:

1. *Spin-up* (periodic). One characteristic event with jump b₀ every T. The
   periodic state just after the event solves y₀ = Φ(T) y₀ + G(T) c + b₀
   (`computePeriodicInitialCondition`). The observed window starts at age
   a ∈ [0, T] of this cycle (default a = T: just before the next
   characteristic event): y(0⁻) = Φ(a) y₀ + G(a) c.
2. *Observed sequence*, clock t ≥ 0: events (t_k, b_k), k = 1..K (b_k from
   `stress_change` columns, as `sys.b`), observation epochs t_j and a
   reference epoch t_ref.

**Data.** Surface displacement at stations, relative to t_ref:

    u(t_j) = H ∫_{t_ref}^{t_j} y dt + u_fixed(t_j) + B o

- H = [KO(:,~locked), LO(:,:,c,1), LO(:,:,c,2)] · P (P = `sys.P`: totals of
  Maxwell + Kelvin rates), stacked over components.
- u_fixed: everything independent of η — the long-term part −H₀ y_pl (t − t_ref),
  the kinematic hinge, coseismic slip of the events (rates jump, displacement
  only through prescribed slip). Supplied by the caller (a helper builds it
  for the example).
- o: nuisance offset per station and component (basis B), eliminated by
  variable projection. With offsets, t_ref drops out (any t_ref ≤ first
  epoch gives the same misfit), so it defaults to the first epoch.

Misfit J = ½ Σ_j r_jᵀ W r_j, r_j = u(t_j) − d_j (after projecting out o).

## 2. Adjoint derivation

Let θ_r = log η_row(r). Since A, c, b all carry the factor 1/η_row(r) on row r,

    ∂(A y + c)/∂θ_r = −e_r ẏ_r(t),   ∂b/∂θ_r = −e_r b_r.

The data enter through ∫_{t_ref}^{t_j} y, so dJ = Σ_j s_jᵀ δ∫_{t_ref}^{t_j} y with
s_j = Hᵀ W r_j: the adjoint source is the **piecewise-constant** density

    s(t) = Σ_{j : t_ref < t < t_j} s_j

(constant between consecutive breakpoints {t_ref, t_j, t_k}). The adjoint λ
solves −λ̇ = Aᵀ λ + s(t) backward from λ(t_end) = 0, continuous through the
events (the jumps are additive and state-independent). Then the gradient over
the observed window is the Stieltjes integral, jumps included:

    ∂J/∂θ_r |_obs = −∫ λ_r(t) dy_r(t)  =  −∫ λ_r ẏ_r dt − Σ_k λ_r(t_k) b_k,r.

**Periodicity.** λ(0) = μ is the sensitivity to the start state. On the spin-up
cycle, λ_p(t) = Φ(a − t)ᵀ μ on [0, a] carries it back to the periodic state,
and the periodicity condition y₀ = Φ(T) y₀ + … gives ν = (I − Φ(T))⁻ᵀ λ_p(0)
and a second weight ψ(t) = Φ(T − t)ᵀ ν on the whole cycle [0, T], whose jump
term is ν_r b₀,r. Altogether

    ∂J/∂θ_r = −[ ∫_obs λ_r dy_r + ∫_0^a λ_p,r ẏ_r dt + ∫_0^T ψ_r ẏ_r dt + ν_r b₀,r ].

**Modal form.** A = V Λ V⁻¹ (one eigendecomposition per η, shared by forward,
adjoint and tangent). With z = V⁻¹ y and adjoint modal coordinates α = Vᵀ λ:

- between events ż_j(t₀+u) = e^{λ_j u} w_j, w = Λ z(t₀) + V⁻¹c — a pure exponential;
- on a segment of length τ with constant source ŝ = Vᵀ s,
  α_i(t₀+u) = e^{λ_i(τ−u)} α_i(τ) + g(λ_i, τ−u) ŝ_i;
- (I − Φ(T))⁻ᵀ is diagonal: 1/(−expm1(λ_i T)), computed with expm1, so slow
  modes keep their digits.

Each time integral is then a pairwise (i, j) closed form,

    ∫_0^τ α_i ż_j du = α_i(τ) w_j E1_ij(τ) + ŝ_i w_j E2_ij(τ),
    E1_ij = ∫_0^τ e^{λ_i(τ−u)} e^{λ_j u} du   (first divided difference of e^{xτ})
    E2_ij = ∫_0^τ g(λ_i, τ−u) e^{λ_j u} du    (second divided difference, nodes λ_i, λ_j, 0)

Summed over all segments (observed and spin-up) and jump terms into one n × n
matrix Q, the gradient is

    ∂J/∂θ = −Re( rowsum( (V⁻ᵀ Q) ∘ V ) )

so the cost is one eigendecomposition + V⁻¹ (already needed by the forward)
+ O(n² · #segments) to fill Q + **one n × n GEMM**. On the example mesh
(Maxwell n = 1061) that is well under the 2 s of the forward; for Burgers
(n = 2103) about 0.1 s.

E1, E2 need a stable evaluation when λ_i ≈ λ_j (and λ → 0): closed forms
with expm1 when |λ_i − λ_j|τ is not small, a Taylor series about the
midpoint below a switch (as `computeModalIntegrals` does for g and h). They
get their own function and a unit test against expm of the 3 × 3 Van Loan
block matrix, including near-equal and complex eigenvalues.

**Checked.** A scratch prototype (random stable 6-state A = D K̃, periodic
spin-up to age a, one observed event, three epochs, t_ref > 0) gives this
gradient equal to central differences of an expm-based forward to 1.3e-9.

### Tangent (forward sensitivity)

For a direction δθ (per row), δy solves δẏ = A δy − δθ ∘ ẏ, with jumps
−δθ ∘ b_k and periodic condition δy₀ = Φ(T) δy₀ + ∫₀ᵀ Φ(T−s)(−δθ∘ẏ) ds − δθ∘b₀.
In modal form the source is B e^{Λu} w with B = V⁻¹ diag(δθ) V, so the same
E1/E2 matrices give δz and δ∫z in closed form. Cost per direction: one n³
product to form B, then O(n²) per segment. For r KL modes the Gauss–Newton
Jacobian H δU costs r GEMMs (r = 10: < 1 s on Maxwell, ~1 s on Burgers).
Tangent and adjoint share every matrix, so the dot-product test
⟨J δθ, δd⟩ = ⟨δθ, Jᵀ δd⟩ checks the implementation, not just the derivation.

## 3. Parameters per cell (Maxwell: 1, Burgers: 3)

The core returns derivatives per **row** of the state; rows map to physical
parameters:

| Rows | Parameter | Gradient |
| --- | --- | --- |
| creeping fault patch | log η′ (rcv.Asigma) | one row each (free by-product) |
| Maxwell, cell m | log η_M,m | sum of its e22 and e23 rows |
| Kelvin, cell m | log η_K,m | sum of its two rows |
| Kelvin, cell m | log μ_K,m | −∫ λ_r (μ_K/η_K) y_r over its two rows (diagonal term, same Q machinery with y instead of ẏ) |

- **Maxwell:** one parameter per cell, log η_M.
- **Burgers:** three per cell. They are not inverted as raw (η_M, η_K, μ_K):
  η_K and μ_K trade off almost perfectly, because postseismic data see the
  Kelvin element mainly through its relaxation time τ_K = η_K/μ_K and its
  strength (amplitude) set by μ_K relative to the elastic modulus G. Default
  parameters per cell:

      p = (log η_M,  log τ_K = log η_K − log μ_K,  log(μ_K/G))

  Because everything is in logs, this is a fixed linear map L of
  (log η_M, log η_K, log μ_K), and the gradient is L⁻ᵀ applied per cell:

      ∂J/∂log η_M = g_M,   ∂J/∂log τ_K = g_ηK,   ∂J/∂log(μ_K/G) = g_ηK + g_μK

  (holding the other two fixed, raising μ_K at fixed τ_K raises η_K too).
  Other choices — e.g. log(η_K/η_M) instead of log τ_K, or one field with
  ratios held fixed — are just a different L; `extractViscosityGradient`
  and `assembleRowPerturbation` take L as an option, with the above as the
  default and 'raw' (log η_M, log η_K, log μ_K) available.
- **Holding one fixed** (user, 2026-11-25): to invert η_M and η_K with μ_K
  known, keep the default parameters and leave `.muKrel` out of the
  directions. With μ_K fixed, log τ_K moves log η_K one-for-one, so the
  ratio parameterisation still holds. Example 6, case B.

## 4. API (three layers)

**Code layout (user, 2026-11-25).** Everything this feature adds lives in a
new top-level folder `inverse/` (sibling of `functions/`), for code that
matters only for inverse or sensitivity work: `solveLinearSequence`,
`computeAdjointGradient`, `computeTangentSequence`, the divided-difference
helper, `extractViscosityGradient`, `assembleRowPerturbation`,
`assembleStationKernels`, `computeDisplacementMisfit`. It calls into
`functions/` (`assembleLinearSystem`, `eigOperator`, `computeModalIntegrals`)
but nothing in `functions/` calls it, and no existing function changes.
Path wiring that comes with the folder:

- `tests/helpers/I2ECPathFixture.m`: add `<root>/inverse`.
- `examples/exampleModel.m` (its `addpath` line): add `inverse/`.
- `docs/conf.py`: add `"inverse"` to `MATLAB_CODE_DIRS`; new
  `docs/api/inverse.rst`; `.github/workflows/docs.yml`: add `inverse/**` to
  the trigger paths.
- When `tasks.md` §7 (`setup.m`) lands on `main`, it should add `inverse/` too.

Naming follows `tasks.md` §4 (verb-first lowerCamelCase), help blocks §9
(H1 line, INPUTS/OUTPUTS, NOTES, SEE ALSO, AUTHORS).

**Forward, cached for the adjoint**

```matlab
seq = struct('T', Trecur, 'b0', sys.b(:,1), 'age', Trecur, ...   % spin-up
             'tEvents', tk, 'bEvents', sys.b(:,2:end), ...          % observed events
             'tObs', tj, 'tRef', tj(1));
fwd = solveLinearSequence(sys, seq);   % fwd.U [n x nObs] = int_tRef^tj y, fwd.Y, and the modes
```

One eigendecomposition (`eigOperator`; balanced or not chosen once by the
same propagation check as `solveAnalyticIvp`, or fixed with `'Balance'`),
V⁻¹, the periodic state in modal form, and all segment states, kept in
`fwd` for the adjoint and tangent. Checked against
`computePeriodicInitialCondition` + `solveAnalyticIvp`.

**Layer 1 — core (no stations)**

```matlab
g   = computeAdjointGradient(fwd, S)      % S [n x nObs] = dJ/dU(:,j); g.eta [n x 1], g.muK [n x 1] per row
dU  = computeTangentSequence(fwd, dTheta) % dTheta [n x r] per-row directions -> dU [n x nObs x r]
gp  = extractViscosityGradient(g, sys)                   % rows -> .fault, .etaM [M x 1], Burgers: .tauK, .muKrel [M x 1]
dTh = assembleRowPerturbation(sys, dp)                   % physical directions -> row directions (same map L)
```

**Layer 2 — observations**

```matlab
obs.H = assembleStationKernels(devl, sys, 'Components', [1 2]);   % [nd x n], = G*P
obs.t = tj; obs.d = d; obs.W = 1./sigma.^2; obs.fixed = ufix; obs.offsets = true;
[J, res] = computeDisplacementMisfit(obs, fwd);   % res.r (projected), res.S = H'*W*r for layer 1
```

**Layer 3 — parameterisation (caller)**: log η = m₀ + Φ ξ; gradient Φᵀ ∇,
Jacobian H · dU(Φ). Example 6 builds Matérn KL modes on cell centroids
inline (no library function; the caller owns this layer).

**Gravity.** The core uses only sys.A, sys.c, sys.b, sys.eta, sys.decay and
sys.map, and the row-scaling rule. Gravity as corrected kernels K_g changes
K̃ only, so nothing changes. If gravity adds state rows with no viscosity
(e.g. accumulated displacement), they get η_row = NaN and are skipped (no
parameter); the adjoint stays valid as long as η enters only as a row scaling.

## 5. Tests (`tests/unit/AdjointTest.m`, small model of `makeSmallModel`)

- Divided differences E1, E2 against 3 × 3 Van Loan expm: separated,
  near-equal (|Δλ|τ = 1e-6 … 1e-1), λ → 0, complex pairs.
- `solveLinearSequence` against `computePeriodicInitialCondition` +
  `solveAnalyticIvp` (periodic state, and the observed window chained by hand).
- Gradient vs central FD (step 1e-3 to 1e-4 in log η, relative error ≲ 1e-6),
  in a random direction and per parameter: **Maxwell**, then **Burgers**
  (default (η_M, τ_K, μ_K/G) and raw parameters).
- Parts separately: a misfit that sees only the start state (periodicity
  adjoint alone), one with events only in the observed window, t_ref > 0.
- Tangent vs FD; dot-product test tangent/adjoint (≤ 1e-10 relative).
- Variable projection: gradient of the offset-projected misfit vs FD.
- Gravity-readiness: a row with η_row = NaN gets no gradient and does not
  break the solve (synthetic extra row).

Slow tier (`ExampleScriptsTest`): example 6 runs and recovers ξ within its
noise level.

## 6. Example 6: `examples/example6_adjointInversion.m`

*As built: Burgers, see "Results" above. The plan below was the first
(Maxwell) version.*

Example mesh, Maxwell. True log η(x) = 1e19 Pa s background with a smooth
low-viscosity wedge anomaly built from 6 Matérn KL modes; one characteristic
earthquake every 500 yr (spin-up), observed window: the event at t = 0 then
20 yr of surface displacement at `m.obs` (e.g. 40 epochs, 1 mm noise, offsets
free). Gauss–Newton (Levenberg–Marquardt damping, small prior on ξ) from the
uniform 1e19 model: misfit per iteration, recovered vs true η in the wedge,
fits at three stations. Also reports gradient/forward cost. Expected run
time: under a minute. A short docs page explains the three layers.

## 7. Open questions (for the user)

1. ~~Burgers/Andrade parameters~~ *Answered 2026-11-25: Maxwell η_M per cell; Burgers three per cell, reparameterised (§3); Andrade skipped for now.*
2. ~~Start of the observed window~~ *Default accepted 2026-11-25: age a of the spin-up cycle, default a = T.*
3. ~~Forward periodic solve~~ *Default accepted: modal in the inversion path; tests compare it with `computePeriodicInitialCondition`.*
4. ~~Fault viscosity~~ *Default accepted: its gradient is returned; example 6 inverts the mantle only.*
5. ~~u_fixed~~ *Default accepted: supplied by the caller (example 6 builds the long-term trend).*
