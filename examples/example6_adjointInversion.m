% EXAMPLE 6: invert surface displacement for a Burgers mantle and fault with adjoint gradients
%
% How to use the adjoint gradients of inverse/ in an inversion. A
% synthetic Burgers mantle (Maxwell dashpot eta_M in series with one
% Kelvin element eta_K, mu_K) and a Burgers fault below the locked zone
% (fault dashpot eta' in series with a Kelvin element: transient
% afterslip) are recovered from GNSS-like displacement time series by
% fmincon, with the gradient of the misfit from one adjoint solve per
% evaluation:
%   - truth: a thermal model of the subduction zone (half-space cooling of
%     the oceanic plate and slab, a hot mantle wedge under a conductive
%     overriding plate) turned into eta_M by an Arrhenius law; viscosity
%     is high in the cold lithosphere, lowest below it and rises again with
%     depth. eta_K = eta_M/10 and one mu_K = G/2 everywhere. The creeping
%     fault (30 to 49 km) is a Burgers fault with eta' = 1e19 Pa s, eta_K =
%     1e18 Pa s and mu_K = G/2, each divided by the fault width as in
%     examples 2-5 (tauK = 2.1 yr);
%   - data: the periodic cycle (one earthquake every 500 yr) as spin-up;
%     21 stations (x and z), every 3 months for the 5 yr before the
%     earthquake and log-spaced for the 10 yr after, 5 mm noise, a few
%     gaps, two stations installed after the earthquake, one that stops;
%     a free offset per station and component. The earthquake source and
%     its timing are known, so its coseismic offsets are fixed, as are the
%     long-term motion and the kinematic hinge;
%   - parameters: in the mantle one uniform value each of log eta_M,
%     log tauK = log(eta_K/mu_K) and log(mu_K/G), plus k Karhunen-Loeve
%     modes (Matern covariance on the mesh) for log eta_M and log tauK,
%     with k = 0, 2, 5, 10, 20; on the fault one uniform value each of
%     log eta', log tauK and log(mu_K/G) (6 to 46 unknowns against ~4000
%     data);
%   - misfit: data only (no regularisation); fmincon with the adjoint
%     gradient and bounds (mantle eta_M in 1e17 to 1e23 Pa s, tauK in 0.01
%     to 1e4 yr; fault eta' in 1e16 to 1e22 Pa s, tauK in 1e-3 to 30 yr;
%     mu_K/G in 0.05 to 20 for both).
% See docs/examples/example6.md and docs/dev/adjoint.md.
%
% Run from any folder:  run examples/example6_adjointInversion
%
% AUTHORS
%   Rishav Mallick, 2026

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; Vpl = m.Vpl; yr = m.yr; G = m.mu;
cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);     % fault, earthquake, loading
locked = cyc.locked; T = cyc.Trecur;

%% the true mantle
th = thermalModel(shz, rcv.dip(1), Vpl, yr);
muKtrue = 0.5;                                        % mu_K/G, one value
truth = [log(th.etaM*1e-6), log(th.etaM/10*1e-6/(muKtrue*G))];   % log eta_M (MPa s), log tauK (s)
fprintf('true mantle: eta_M = 10^%.1f to 10^%.1f Pa s, tauK = %.2g to %.2g yr, mu_K = G/2\n', ...
    log10(min(th.etaM)), log10(max(th.etaM)), min(exp(truth(:,2)))/yr, max(exp(truth(:,2)))/yr);

%% the true fault: Burgers below the locked zone
% bulk values (Pa s, MPa) divided by the fault width h, as in examples 2-5:
% eta' and eta_K in MPa s/m, mu_K in MPa/m; the fault parameters are
% [log eta', log tauK (s), log(mu_K/k0)] with k0 = G/h
h = sum(rcv.W);
faultTrue = [log(1e-6*1e19/h), log(1e-6*1e18/(0.5*G)), log(0.5)];
fprintf('true fault: eta'' = %.0e Pa s, eta_K = %.0e Pa s, mu_K = G/2 (tauK = %.2f yr), per fault width\n', ...
    exp(faultTrue(1))*h*1e6, exp(faultTrue(2))*0.5*G*1e6, exp(faultTrue(2))/yr);

%% observations: 5 yr before the earthquake, 10 yr after
tEq = 5*yr;                                           % earthquake, observed clock
tObs = [(0:0.25:4.75)*yr, tEq + logspace(log10(1/365), 1, 80)*yr];
sta = 1:20:size(m.obs, 1);                            % 21 stations, every 30 km
xsta = m.obs(sta, 1);
ns = numel(sta); comp = [1 2];
sigma = 5e-3;                                        % m
% missing data: stations 15 and 17 installed 0.5 yr after the earthquake,
% station 3 stops 7 yr after it, station 9 has a 2-yr gap, 3% dropouts
rng(2);
miss = rand(ns, numel(tObs)) < 0.03;
miss([15 17], :) = miss([15 17], :) | tObs < tEq + 0.5*yr;
miss(3, :) = miss(3, :) | tObs > tEq + 7*yr;
miss(9, :) = miss(9, :) | (tObs > tEq + 2*yr & tObs < tEq + 4*yr);
miss = [miss; miss];                                  % both components

% displacement that does not depend on the mantle (all known): the
% long-term motion (plate rate on the fault, corner flow in the mantle) and
% the kinematic hinge since the first epoch, as in example 2, and the
% coseismic slip of the locked zone
xpl = [rcv.Vpl; shz.e22pl; shz.e23pl];
Ghinge = cat(3, m.Ghx(sta,:), m.Ghz(sta,:));
fixed = [];
for c = comp
    Gc = [m.devl.KO(sta,:,c), m.devl.LO(sta,:,c,1), m.devl.LO(sta,:,c,2)];
    vlong = -Gc*xpl - Ghinge(:,:,c)*(m.hinge.Vpl*Vpl);
    ucos = m.devl.KO(sta, locked, c) * cyc.slip(locked);
    fixed = [fixed; vlong*(tObs - tObs(1)) + ucos*(tObs >= tEq)]; %#ok<AGROW>
end

% the time line: the observed window starts 5 yr before the end of the
% periodic cycle, and the earthquake at tEq is that cycle's next one
ctx = struct('G', G, 'h', h, 'evl', m.evl, 'rcv', rcv, 'shz', shz, 'locked', locked, ...
    'sc', cyc.stress_change, 'T', T, 'age', T - tEq, 'tEq', tEq, 'tObs', tObs);

%% synthetic data
model = struct('cells', [truth, log(muKtrue)*ones(shz.N, 1)], 'fault', faultTrue);
fwdTrue = solveFields(model, ctx);
afterslipTrue = meanAfterslip(fwdTrue, ctx);
obs.H = assembleStationKernels(m.devl, fwdTrue.sys, 'Components', comp, 'Stations', sta);
obs.W = 1/sigma^2; obs.fixed = fixed; obs.offsets = true;
rng(1);
obs.d = obs.H*fwdTrue.U + fixed + 0.1*randn(2*ns, 1) + sigma*randn(2*ns, numel(tObs));
obs.d(miss) = NaN;
nd = nnz(~miss);
fprintf('%d stations x 2 components, %d epochs: %d data (%d missing); %d Burgers states\n', ...
    ns, numel(tObs), nd, nnz(miss), numel(fwdTrue.sys.eta));

%% parameters: uniform values plus k KL modes for log eta_M and log tauK
w = shz.area / sum(shz.area);
d = sqrt((shz.xc(:,1) - shz.xc(:,1)').^2 + (shz.xc(:,2) - shz.xc(:,2)').^2) / 60e3;
C = (1 + sqrt(3)*d) .* exp(-sqrt(3)*d);               % Matern (nu = 3/2), 60 km
[U, L] = eig(sqrt(w) .* C .* sqrt(w'), 'vector');
[~, order] = sort(L, 'descend');
KL = U(:, order(1:20)) ./ sqrt(w);                    % the 20 smoothest modes, unit rms
% start: mantle eta_M = 1e19 Pa s, eta_K = 1e18 Pa s, mu_K = G (tauK ~ 1 yr);
% fault eta' = 1e18 Pa s, eta_K = 1e18 Pa s, mu_K = G (tauK ~ 1 yr)
start = [log(1e19*1e-6), log(1e18*1e-6/G), 0, log(1e-6*1e18/h), log(1e18*1e-6/G), 0];
% bounds: mantle eta_M (MPa s), tauK (s), mu_K/G; fault eta' (MPa s/m),
% tauK (s), mu_K/k0
bounds = [log([1e17 1e23]*1e-6); log([0.01 1e4]*yr); log([0.05 20]); ...
          log([1e16 1e22]*1e-6/h); log([1e-3 30]*yr); log([0.05 20])];

%% inversions for k = 0, 2, 5, 10, 20
kList = [0 2 5 10 20];
res = cell(size(kList));
p = start(:);                                         % [cM; cT; cU; fault] for k = 0
for i = 1:numel(kList)
    k = kList(i);
    par = parameterisation(KL(:, 1:k), bounds);
    if i > 1                                          % warm start: new modes at zero
        p = par.from(res{i-1}.p, kList(i-1));
    end
    if i == 1
        checkGradient(p, par, obs, ctx, nd);
    end
    res{i} = invert(p, par, obs, ctx, nd);
    p = res{i}.p;
    q = par.model(p);
    fprintf(['  k = %d (%2d unknowns): chi^2/N = %7.3f after %3d evaluations (%3.0f s); ' ...
        'mantle mu_K/G = %.2f; fault eta'' = %.2g Pa s, tauK = %.3f yr, mu_K/G = %.2f\n'], ...
        k, numel(p), res{i}.chi, res{i}.nfun, res{i}.time, exp(q.cells(1,3)), ...
        exp(q.fault(1))*h*1e6, exp(q.fault(2))/yr, exp(q.fault(3)));
end
fprintf('  true: mantle mu_K/G = %.2f; fault eta'' = %.2g Pa s, tauK = %.3f yr, mu_K/G = %.2f\n', ...
    muKtrue, exp(faultTrue(1))*h*1e6, exp(faultTrue(2))/yr, exp(faultTrue(3)));

%% figures
st = exampleStyle();
shift = [6, -log10(yr)];                              % MPa s -> Pa s, s -> yr
rowName = {'log_{10}(\eta_M / Pa s)', 'log_{10}(\tau_K / yr)'};
st.figure(1, 1500, 420)
tl = tiledlayout(2, numel(kList) + 1, 'TileSpacing', 'tight', 'Padding', 'compact');
for f = 1:2
    lim = prctile(truth(:,f)/log(10) + shift(f), [1 90]);
    for i = 0:numel(kList)
        ax = nexttile;
        if i == 0
            v = truth(:,f); ttl = 'true';
        else
            fi = parameterisation(KL(:, 1:kList(i)), bounds).model(res{i}.p);
            v = fi.cells(:,f); ttl = sprintf('k = %d, \\chi^2/N = %.2f', kList(i), res{i}.chi);
        end
        meshPlot(shz, rcv, v/log(10) + shift(f), xsta, st)
        clim(lim), colormap(ax, turbo(256))
        if f == 1, title(ttl), end
        if i == 0, ylabel(rowName{f}), else, set(ax, 'YTickLabel', []), end
        if f < 2, set(ax, 'XTickLabel', []), end
        if i == numel(kList), colorbar; end
    end
end
xlabel(tl, 'distance from the trench (km); depth (km) on the vertical axes; triangles: stations', ...
    'FontSize', 11, 'Color', st.ink2)
st.axes()

% data and fits less the known terms (coseismic step, long-term motion,
% hinge) and the fitted offset, i.e. the mantle and afterslip signal, at
% five stations: data as filled circles, the first and last k as dashed and
% solid lines. Interseismic (5 yr before the earthquake) and postseismic
% (10 yr after, log time) separately, each shifted to start at 0; the
% interseismic series are spread apart by a constant per station so the
% dots and fits do not overlap
show = [3 5 8 13 21];                                 % stations at x = -40, 20, 110, 260, 500 km
colours = [st.series; 0.55 0.36 0.85; 0.80 0.60 0.10];
before = tObs < tEq; after = ~before;
i0 = find(before, 1, 'last');                         % last epoch before the earthquake
windows = {before, after}; names = {'Interseismic', 'Postseismic'};
for f = 1:2
    st.figure(f + 1, ifelse(f == 1, 1100, 1300), 420)
    tiledlayout(1, ifelse(f == 1, 2, 3), 'TileSpacing', 'compact', 'Padding', 'compact');
    sel = windows{f};
    t = (tObs(sel) - tEq)/yr;
    for c = 1:2
        nexttile, hold on
        for j = 1:numel(show)
            row = show(j) + (c - 1)*ns;
            col = colours(j,:);
            % each model less its own fitted offset, shifted to 0 at the
            % first epoch (interseismic) or the last before the earthquake
            % (postseismic); the data with the last model's offset and shift
            iref = ifelse(f == 1, 1, i0);
            signal = @(i) res{i}.u(row,:) - obs.fixed(row,:) - res{i}.offset(row);
            ref = signal(numel(kList)); ref = ref(iref);
            spread = ifelse(f == 1, (j - 1)*5e-2, 0);   % m
            dm = obs.d(row,:) - obs.fixed(row,:) - res{end}.offset(row) - ref + spread;
            plot(t, 1e2*dm(sel), 'o', 'MarkerSize', 3, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
                'DisplayName', sprintf('x = %.0f km', xsta(show(j))/1e3))
            for i = [1 numel(kList)]
                ui = signal(i); ui = ui - ui(iref) + spread;
                plot(t, 1e2*ui(sel), ifelse(i == 1, '--', '-'), ...
                    'Color', col, 'LineWidth', ifelse(i == 1, 1, 1.8), 'HandleVisibility', 'off')
            end
        end
        % style key in grey: data, the first and the last model
        plot(NaN, NaN, 'o', 'MarkerSize', 3, 'MarkerFaceColor', st.ink2, 'MarkerEdgeColor', 'none', 'DisplayName', 'data')
        plot(NaN, NaN, '--', 'Color', st.ink2, 'LineWidth', 1, 'DisplayName', sprintf('model, k = %d', kList(1)))
        plot(NaN, NaN, '-', 'Color', st.ink2, 'LineWidth', 1.8, 'DisplayName', sprintf('model, k = %d', kList(end)))
        grid on
        if f == 1, xlim([-5 0]), else, set(gca, 'XScale', 'log'), xlim([1/365 10]), end
        if f == 1 && c == 2
            legend('Location', 'eastoutside', 'Box', 'off', 'FontSize', 8)
        elseif f == 2 && c == 1
            legend('Location', 'best', 'Box', 'off', 'FontSize', 8, 'NumColumns', 2)
        end
        xlabel('time since the earthquake (yr)')
        ylabel(sprintf('%s displacement (cm)%s', ifelse(c == 1, 'horizontal', 'vertical'), ...
            ifelse(f == 1, ', offset per station', '')))
        title(sprintf('%s, %s', names{f}, ifelse(c == 1, 'horizontal', 'vertical')))
    end
    if f == 2
        nexttile, hold on
        plot(t, 1e2*afterslipTrue(sel), 'Color', st.ink, 'LineWidth', 2, 'DisplayName', 'true')
        for i = [1 numel(kList)]
            plot(t, 1e2*res{i}.afterslip(sel), ifelse(i == 1, '--', '-'), ...
                'Color', st.series(ifelse(i == 1, 2, 1), :), 'LineWidth', 1.5, 'DisplayName', sprintf('k = %d', kList(i)))
        end
        set(gca, 'XScale', 'log'), grid on, xlim([1/365 10])
        legend('Location', 'northwest', 'Box', 'off', 'FontSize', 8)
        xlabel('time since the earthquake (yr)'), ylabel('mean afterslip (cm)')
        title('Afterslip below the locked zone')
    end
    st.axes()
end

st.figure(4, 760, 300)
meshPlot(shz, rcv, th.T, xsta, st), colormap(gca, turbo(256))
cb = colorbar; cb.Label.String = 'temperature (C)';
xlabel('distance from the trench (km)'), ylabel('depth (km)')
title('Thermal model of the true mantle')
st.axes()

%% ---------------------------------------------------------------------
function par = parameterisation(Phi, bounds)
% p = [cM; aM (k); cT; aT (k); cU; fE; fT; fU]:
%   mantle  log eta_M = cM + Phi*aM,  log tauK = cT + Phi*aT,  log(mu_K/G) = cU;
%   fault   log eta' = fE,  log tauK = fT,  log(mu_K/k0) = fU (uniform).
% Returns the map p -> model (cell fields and fault values), the chain
% rule for the gradient (uniform fault values: sums over the patches), and
% the bounds: linear inequalities on the mantle fields, boxes on the rest.
[M, k] = size(Phi);
B = [ones(M, 1), Phi];                                % one field from (c, a)
iM = 1:k+1; iT = k+1 + (1:k+1); iU = 2*k + 3; iF = 2*k + 3 + (1:3);
n = 2*k + 6;
par.k = k;
par.model = @(p) struct('cells', [B*p(iM), B*p(iT), p(iU)*ones(M, 1)], 'fault', p(iF)');
par.chain = @(gp) [B'*gp.etaM; B'*gp.tauK; sum(gp.muKrel); ...
                   sum(gp.fault); sum(gp.faultTauK); sum(gp.faultMuKrel)];
Z = zeros(M, k+1); z = zeros(M, 4);
par.A = [B, Z, z; -B, Z, z; Z, B, z; Z, -B, z];
par.b = [bounds(1,2)*ones(M, 1); -bounds(1,1)*ones(M, 1); ...
         bounds(2,2)*ones(M, 1); -bounds(2,1)*ones(M, 1)];
par.lb = -Inf(n, 1); par.ub = Inf(n, 1);
par.lb([iU iF]) = bounds(3:6,1); par.ub([iU iF]) = bounds(3:6,2);
% coefficients of a smaller k in this layout, new modes at zero
par.from = @(q, kq) [q(1:kq+1); zeros(k - kq, 1); q(kq+1 + (1:kq+1)); zeros(k - kq, 1); q(end-3:end)];
end

function [f, g] = objective(p, par, obs, ctx, nd)
% chi^2/N of the data and its gradient: one forward and one adjoint solve
try
    fwd = solveFields(par.model(p), ctx);
catch err                                              % a mode that no longer decays
    if ~startsWith(err.identifier, 'I2EC:'), rethrow(err), end
    f = Inf; g = zeros(size(p));
    return
end
[J, r] = computeDisplacementMisfit(obs, fwd);
f = 2*J/nd;
if nargout > 1
    gp = extractViscosityGradient(computeAdjointGradient(fwd, r.S), fwd.sys);
    g = 2*par.chain(gp)/nd;
end
end

function out = invert(p, par, obs, ctx, nd)
% fmincon with the adjoint gradient and the bounds on the fields
opts = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, 'MaxIterations', 40, ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8, 'Display', 'off');
t0 = tic;
[p, ~, ~, info] = fmincon(@(q) objective(q, par, obs, ctx, nd), p, par.A, par.b, ...
    [], [], par.lb, par.ub, [], opts);
fwd = solveFields(par.model(p), ctx);
[J, r] = computeDisplacementMisfit(obs, fwd);
out = struct('p', p, 'chi', 2*J/nd, 'u', r.u, 'offset', r.offset, ...
    'afterslip', meanAfterslip(fwd, ctx), 'nfun', info.funcCount, 'time', toc(t0));
end

function checkGradient(p, par, obs, ctx, nd)
% adjoint gradient against a central difference along a random direction
rng(3);
dp = randn(size(p)); h = 1e-4;
[~, g] = objective(p, par, obs, ctx, nd);
fd = (objective(p + h*dp, par, obs, ctx, nd) - objective(p - h*dp, par, obs, ctx, nd))/(2*h);
fprintf('gradient check: adjoint %.6e, central difference %.6e (relative difference %.1e)\n', ...
    g'*dp, fd, abs(g'*dp - fd)/abs(fd));
end

function fwd = solveFields(q, ctx)
% forward model for q.cells = [log eta_M (MPa s), log tauK (s),
% log(mu_K/G)] per cell and q.fault = [log eta' (MPa s/m), log tauK (s),
% log(mu_K/k0)], k0 = G/h: Burgers mantle and fault; the periodic cycle,
% then the observed window from 5 yr before its next earthquake
shz = ctx.shz; rcv = ctx.rcv;
shz.alpha = exp(-q.cells(:,1));                       % 1/eta_M
shz.muK = ctx.G * exp(q.cells(:,3));                  % mu_K = G * (mu_K/G)
shz.etaK = exp(q.cells(:,2)) .* shz.muK;              % eta_K = tauK * mu_K
rcv.Asigma = exp(q.fault(1)) * ones(rcv.N, 1);        % eta'
rcv.muK = ctx.G/ctx.h * exp(q.fault(3)) * ones(rcv.N, 1);
rcv.etaK = exp(q.fault(2)) * rcv.muK;                 % eta_K = tauK * mu_K
sys = assembleLinearSystem(ctx.evl, rcv, shz, ctx.locked, ctx.sc);
seq = struct('T', ctx.T, 'b0', sys.b, 'age', ctx.age, 'tEvents', ctx.tEq, ...
    'bEvents', sys.b, 'tObs', ctx.tObs, 'tRef', ctx.tObs(1));
fwd = solveLinearSequence(sys, seq, 'Balance', true);  % balanced: viscosity spans decades
end

function th = thermalModel(shz, dipDeg, Vpl, yr)
% Temperature and Maxwell viscosity of the mesh from 2-D analytic thermal
% solutions (Turcotte & Schubert, Geodynamics, ch. 4 and 6): half-space
% cooling of a 20 Myr plate below the seafloor and the slab top; in the
% mantle wedge a hot core with a boundary layer on the slab,
% erf(d_s/sqrt(kappa*s/V)), under a conductive overriding plate,
% erf(z/40 km); an adiabat of 0.4 K/km. Viscosity: diffusion creep,
% E = 300 kJ/mol, Vact = 12 cm^3/mol, 1e19 Pa s at 100 km depth, capped at
% 1e22 Pa s: high in the lithosphere, low below, rising again with depth.
kappa = 1e-6; Tp = 1350; age = 20e6*yr; Lo = 40e3; dTdz = 0.4e-3;
E = 300e3; Vact = 12e-6; R = 8.314; rho = 3300; g = 9.8;
x = shz.xc(:,1); z = -shz.xc(:,2);                   % depth positive
dip = dipDeg*pi/180;
above = x*sin(dip) - z*cos(dip) > 0;                  % mantle wedge
dSlab = abs(x*sin(dip) - z*cos(dip));                 % distance to the slab top
s = max(x*cos(dip) + z*sin(dip), 1e3);                % down-dip distance from the trench
dOc = dSlab; dOc(x < 0) = z(x < 0);                   % below the seafloor or the slab top
T = Tp * erf(dOc / (2*sqrt(kappa*age)));
T(above) = Tp * erf(dSlab(above) ./ sqrt(kappa*s(above)/Vpl)) .* erf(z(above)/Lo);
T = T + dTdz*z;
Tr = Tp + dTdz*100e3; Pr = rho*g*100e3;
th.T = T;
th.etaM = 1e19 * exp((E + rho*g*z*Vact)./(R*(T + 273)) - (E + Pr*Vact)/(R*(Tr + 273)));
th.etaM = min(th.etaM, 1e22);                         % Pa s
end

function s = meanAfterslip(fwd, ctx)
% width-weighted mean slip of the creeping patches since the first epoch,
% less creep at the plate rate (dashpot plus Kelvin elements)
slip = extractRatesFromLinearState(fwd.U, fwd.sys);
c = ~ctx.locked; W = ctx.rcv.W(c);
s = W'*slip(c,:)/sum(W) - ctx.rcv.Vpl(1)*(ctx.tObs - ctx.tObs(1));
end

function meshPlot(shz, rcv, v, xsta, st)
% mantle mesh coloured by v, the fault and the stations, undistorted
plotShz2d(shz, v), hold on, plotPatches2d(rcv)
plot(xsta/1e3, 6*ones(size(xsta)), 'v', 'Color', st.ink, 'MarkerFaceColor', st.ink, 'MarkerSize', 4)
axis equal, xlim([-100 500]), ylim([-200 12])
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
