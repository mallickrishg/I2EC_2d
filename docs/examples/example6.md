# Example 6: Burgers mantle and fault inversion with adjoint gradients

How to use the gradients of `inverse/` in an inversion. A synthetic Burgers
mantle and a Burgers fault are recovered from GNSS-like displacement time
series with MATLAB's `fmincon`, and each evaluation of the misfit and its
gradient costs one forward and one adjoint solve. The example shows the
workflow and how the fit improves as the parameterisation gets richer; it is
not a study of resolution. Run it with

```matlab
run examples/example6_adjointInversion
```

It takes about 16 minutes: five inversions of 50 to 110 evaluations each, at
about 2 s per evaluation (2122 states).

## 1. The synthetic truth

**Mantle.** A thermal model of the subduction zone (half-space cooling of a
20 Myr plate and of the slab, a hot wedge below a conductive overriding plate,
an adiabat) gives $\eta_M$ through an Arrhenius law: high in the cold
lithosphere and slab, lowest below them, rising again with depth
($10^{18.9}$ to $10^{22}$ Pa s). Each cell is a Burgers body with
$\eta_K = \eta_M/10$ and one $\mu_K = G/2$ everywhere.

```{figure} ../_static/examples/example6_adjointInversion_4.png
:alt: Temperature of the true mantle on the mesh
:width: 70%

Temperature of the true mantle, from which $\eta_M$ follows.
```

**Fault.** The fault is locked above 30 km and creeps below it (30 to 49 km).
The creeping section is a Burgers fault: the fault dashpot $\eta'$ in series
with a Kelvin element (transient afterslip), with $\eta' = 10^{19}$ Pa s,
$\eta_K = 10^{18}$ Pa s and $\mu_K = G/2$, each divided by the fault width as
in examples 2–5 ($\tau_K = \eta_K/\mu_K = 2.1$ yr). See
[the fault](../rheology.md) for the units.

**Earthquake and loading.** One earthquake every 500 years releases the slip
deficit of the locked zone. Its source and timing are known, so its
coseismic displacement is a fixed term, as are the long-term motion and the
kinematic hinge.

## 2. The data

21 stations every 30 km, both components, 5 years before the earthquake
every 3 months and 10 years after it at 80 log-spaced epochs, with 5 mm noise
and a free offset per station and component. The record is messy on purpose:
two stations are installed 6 months after the earthquake, one stops after 7
years, one has a 2-year gap, and 3% of the epochs are missing at random
(3788 data, 412 missing). Missing data are `NaN` in `obs.d`; they carry no
misfit and no adjoint source.

## 3. The forward model and the misfit

{mat:func}`solveLinearSequence` solves the periodic cycle directly (no
spin-up cycles) and then the observed window, which starts 5 years before the
cycle's next earthquake:

```matlab
sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change);
seq = struct('T', T, 'b0', sys.b, 'age', T - tEq, ...   % periodic cycle
             'tEvents', tEq, 'bEvents', sys.b, ...       % observed earthquake
             'tObs', tObs, 'tRef', tObs(1));             % epochs, reference
fwd = solveLinearSequence(sys, seq);

obs.H = assembleStationKernels(devl, sys, 'Stations', sta);
obs.d = d; obs.W = 1/sigma^2; obs.fixed = fixed; obs.offsets = true;
[J, res] = computeDisplacementMisfit(obs, fwd);
```

`obs.fixed` holds everything that does not depend on the viscosities, here
the coseismic step, the long-term motion and the hinge, so they cancel
exactly in the residual. The misfit is the data alone: no regularisation.

## 4. Parameters and the gradient

The data see a Kelvin element mainly through its relaxation time and its
strength, so each Burgers body is parameterised by logs of ratios,

$$
\log\eta_M \;(\text{or } \log\eta'), \qquad \log\tau_K = \log(\eta_K/\mu_K),
\qquad \log(\mu_K/G),
$$

and the parameterisation is the caller's (built in the script):

| Unknowns | Mantle | Fault (creeping section) |
| --- | --- | --- |
| $\log\eta_M$ or $\log\eta'$ | uniform value + $k$ KL modes | uniform |
| $\log\tau_K$ | uniform value + $k$ KL modes | uniform |
| $\log(\mu_K/G)$ | uniform | uniform |

The KL modes are the smoothest eigenvectors of a Matérn ($\nu = 3/2$, 60 km)
covariance on the mesh, for $k = 0, 2, 5, 10, 20$: 6 to 46 unknowns. The
gradient goes from state rows to cells and patches
({mat:func}`extractViscosityGradient`) and then to the unknowns by the chain
rule; a uniform fault value is the sum over the creeping patches:

```matlab
g  = computeAdjointGradient(fwd, res.S);          % per state row
gp = extractViscosityGradient(g, sys);            % per cell and patch
grad = [B'*gp.etaM; B'*gp.tauK; sum(gp.muKrel); ...          % B = [1, KL modes]
        sum(gp.fault); sum(gp.faultTauK); sum(gp.faultMuKrel)];
```

`fmincon` (interior point) takes this gradient, with bounds on every field:
$\eta_M$ in $10^{17}$ to $10^{23}$ Pa s and $\tau_K$ in 0.01 to $10^4$ yr in
the mantle, $\eta'$ in $10^{16}$ to $10^{22}$ Pa s and $\tau_K$ in
$10^{-3}$ to 30 yr on the fault, $\mu_K/G$ in 0.05 to 20 for both. Each $k$
starts from the previous solution with the new modes at zero. The script
first checks the gradient against a central difference:

```text
gradient check: adjoint -7.532525e+01, central difference -7.532525e+01 (relative difference 1.9e-08)
```

## 5. Results

```text
  k = 0 ( 6 unknowns): chi^2/N = 206.270 after  54 evaluations (114 s); mantle mu_K/G = 1.46; fault eta' = 2.2e+20 Pa s, tauK = 1.868 yr, mu_K/G = 0.52
  k = 2 (10 unknowns): chi^2/N =  43.462 after  72 evaluations (158 s); mantle mu_K/G = 0.20; fault eta' = 1.6e+19 Pa s, tauK = 1.984 yr, mu_K/G = 0.51
  k = 5 (16 unknowns): chi^2/N =  12.011 after  95 evaluations (205 s); mantle mu_K/G = 0.58; fault eta' = 8.3e+18 Pa s, tauK = 2.238 yr, mu_K/G = 0.48
  k = 10 (26 unknowns): chi^2/N =   4.983 after  98 evaluations (218 s); mantle mu_K/G = 0.49; fault eta' = 8.4e+18 Pa s, tauK = 2.203 yr, mu_K/G = 0.49
  k = 20 (46 unknowns): chi^2/N =   2.776 after 110 evaluations (241 s); mantle mu_K/G = 0.77; fault eta' = 1.1e+19 Pa s, tauK = 2.306 yr, mu_K/G = 0.46
  true: mantle mu_K/G = 0.50; fault eta' = 1e+19 Pa s, tauK = 2.116 yr, mu_K/G = 0.50
```

The fit improves steadily with $k$, from $\chi^2/N = 206$ for uniform fields
to 2.8 with 20 modes. The fault is recovered from $k = 5$ on, within about
20% in $\eta'$ and 10% in $\tau_K$ and $\mu_K$, and its afterslip from
$k = 0$. The mantle fields that fit the data are not unique: they change
shape with $k$ and do not reproduce the thin cold layers of the truth, which
surface data barely constrain. Only the broad pattern, low viscosity and
short $\tau_K$ below the fore-arc at 80 to 150 km depth, persists.

```{figure} ../_static/examples/example6_adjointInversion_1.png
:alt: True and recovered log eta_M and log tauK on the mesh for each k
:width: 100%

$\log_{10}\eta_M$ (top) and $\log_{10}\tau_K$ (bottom): the truth (left) and
the fields recovered with $k = 0, 2, 5, 10, 20$ modes. Triangles are the
stations.
```

```{figure} ../_static/examples/example6_adjointInversion_2.png
:alt: Interseismic displacement at five stations, data and fits
:width: 80%

The 5 years before the earthquake at stations 40 km seaward of the trench and
20, 110, 260 and 500 km landward, with the fixed terms and fitted offsets
removed and 5 cm between stations: data (dots), uniform fields ($k = 0$,
dashed) and $k = 20$ (solid).
```

```{figure} ../_static/examples/example6_adjointInversion_3.png
:alt: Postseismic displacement at five stations and the afterslip, data and fits
:width: 100%

The 10 years after the earthquake (log time) at the same stations, and the
mean afterslip of the creeping section, true and recovered. The early
transient is mostly the fault's Kelvin element; the later years separate
$k = 0$ from $k = 20$.
```

## How the gradient is computed

Every row of the operator, the forcing and the coseismic jump is divided by
the viscosity of its element (fault dashpot, Maxwell or Kelvin dashpot), so
the derivative with respect to $\log\eta_r$ is $-\int \lambda_r\,\mathrm{d}y_r$:
the adjoint $\lambda$ against the change of the state, jumps included. The
Kelvin spring enters only the diagonal, so the derivative with respect to
$\log\mu_K$ is $-(\mu_K/\eta_K)\int\lambda_r y_r\,\mathrm{d}t$. The periodic
spin-up adds the adjoint of the periodicity condition, which is diagonal in
modal form. All of these are closed-form in the modes of $A$, which the
forward has already computed, so the gradient costs about one more forward
solve, whatever the number of unknowns or the spacing of the epochs. The
tangent ({mat:func}`computeTangentSequence`, one solve per direction) gives
Jacobian columns for Gauss–Newton methods. The derivation, the tests and the
cost are in `docs/dev/adjoint.md`.
