% EXAMPLE 6: viscosity inversion with adjoint gradients and tangent Gauss-Newton
%
% A synthetic inversion for the mantle viscosity eta(x) of the example
% model (Maxwell, as example 2) from surface-displacement time series:
%   - forward model: the periodic cycle (one earthquake every 500 yr) as
%     spin-up, then an observed earthquake at t = 0 and 20 yr of surface
%     displacement at 21 stations (x and z), with 5 mm noise and an unknown
%     offset per station and component;
%   - parameters: log eta = log(1e19 Pa s) + Phi*xi, with Phi the first 6
%     Karhunen-Loeve modes of a Matern covariance on the mantle mesh
%     (the caller's parameterisation layer, built here);
%   - Gauss-Newton (Levenberg-Marquardt) on xi with a unit Gaussian prior;
%     the Jacobian comes from 6 tangent solves that reuse the forward's
%     eigendecomposition, and the adjoint gradient is checked against it.
% The inverse/ functions work per mantle cell and state row; only the last
% layer (Phi) is specific to this script. See docs/dev/adjoint.md.
%
% Run from any folder:  run examples/example6_adjointInversion
%
% AUTHORS
%   Rishav Mallick, 2026

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; Vpl = m.Vpl; yr = m.yr;
cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);     % rheology, event, loading
locked = cyc.locked; T = cyc.Trecur;
eta0 = 1e19*1e-6;                                    % MPa s

%% parameterisation: Matern (nu = 3/2) KL modes on the mesh
ell = 60e3; r = 6;
d = sqrt((shz.xc(:,1) - shz.xc(:,1)').^2 + (shz.xc(:,2) - shz.xc(:,2)').^2);
C = (1 + sqrt(3)*d/ell) .* exp(-sqrt(3)*d/ell);      % unit variance in log eta
w = shz.area / sum(shz.area);
[U, L] = eig(sqrt(w) .* C .* sqrt(w'), 'vector');
[L, order] = sort(L, 'descend');
Phi = (U(:, order(1:r)) ./ sqrt(w)) .* sqrt(L(1:r))';   % [M x r], KL modes
fprintf('%d KL modes hold %.0f%% of the prior variance\n', r, 100*sum(L(1:r))/sum(L));
xiTrue = [-1.2; 0.8; 0.6; -0.5; 0.4; -0.3];

%% observations: 21 stations, x and z, 40 epochs over 20 yr
sta = 1:20:size(m.obs, 1);
tObs = linspace(0.5, 20, 40)*yr;
comp = [1 2];
sigma = 5e-3;                                        % m
% displacement that does not depend on viscosity: the long-term state
% (plate rate on the fault, corner flow in the mantle) and the kinematic
% hinge, accumulated from the first epoch (see example 2)
xpl = [rcv.Vpl; shz.e22pl; shz.e23pl];
Ghinge = cat(3, m.Ghx(sta,:), m.Ghz(sta,:));
vfix = [];
for c = comp
    Gc = [m.devl.KO(sta,:,c), m.devl.LO(sta,:,c,1), m.devl.LO(sta,:,c,2)];
    vfix = [vfix; -Gc*xpl - Ghinge(:,:,c)*(m.hinge.Vpl*Vpl)]; %#ok<AGROW>
end
fixed = vfix * (tObs - tObs(1));

% forward model for a given xi: system, time line, solve
forward = @(xi, bal) solveAt(xi, bal, eta0, Phi, m.evl, rcv, shz, locked, cyc.stress_change, T, tObs);

%% synthetic data
[fwdTrue, balance] = forward(xiTrue, 'auto');
obs.H = assembleStationKernels(m.devl, fwdTrue.sys, 'Components', comp, 'Stations', sta);
obs.W = 1/sigma^2; obs.fixed = fixed; obs.offsets = true;
rng(1);
obs.d = obs.H*fwdTrue.U + fixed + 0.1*randn(size(obs.H, 1), 1) + sigma*randn(size(fixed));
nd = numel(obs.d);

%% Gauss-Newton (Levenberg-Marquardt) from the uniform model
xi = zeros(r, 1); damp = 1e-2; maxit = 10;
hist = zeros(0, 2);               % data misfit chi^2/N, model error
for it = 0:maxit
    tic
    fwd = forward(xi, balance);
    dTheta = assembleRowPerturbation(fwd.sys, struct('etaM', Phi));   % KL modes -> state rows
    dU = computeTangentSequence(fwd, dTheta);                          % r tangent solves
    [J, res] = computeDisplacementMisfit(obs, fwd, dU);
    hist(it+1, :) = [2*J/nd, norm(xi - xiTrue)];
    if it == 0
        % adjoint gradient (rows -> cells -> KL modes) against the Jacobian
        g = computeAdjointGradient(fwd, res.S);
        gp = extractViscosityGradient(g, fwd.sys);
        gXi = Phi' * gp.etaM;
        fprintf('adjoint vs tangent gradient: relative difference %.1e\n', ...
            norm(gXi - res.Jw'*res.rw)/norm(gXi));
        J0 = res;
    end
    fprintf('iteration %d: chi^2/N = %7.3f, |xi - xi_true| = %.3f  (%.1f s)\n', ...
        it, hist(it+1, 1), hist(it+1, 2), toc);
    % LM step on 1/2|rw + Jw*dxi|^2 + 1/2|xi + dxi|^2, damped until it improves;
    % stop when the objective no longer decreases by 1e-6
    obj = J + 0.5*(xi'*xi);
    improved = false;
    while damp < 1e6 && it < maxit
        dxi = -[res.Jw; eye(r); sqrt(damp)*eye(r)] \ [res.rw; xi; zeros(r, 1)];
        Jt = computeDisplacementMisfit(obs, forward(xi + dxi, balance));
        objt = Jt + 0.5*norm(xi + dxi)^2;
        if objt < obj
            improved = obj - objt > 1e-6*obj;
            break
        end
        damp = 10*damp;
    end
    if ~improved, break, end
    xi = xi + dxi; damp = max(damp/10, 1e-6);
end
nit = size(hist, 1) - 1;
fprintf('xi true:      %s\nxi recovered: %s\n', mat2str(xiTrue', 2), mat2str(xi', 2));

%% figures
st = exampleStyle();
logTrue = Phi*xiTrue/log(10);                     % log10(eta / 1e19 Pa s)
logRec = Phi*xi/log(10);
st.figure(1, 1000, 330)
tl = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
for k = 1:2
    nexttile
    plotShz2d(shz, ifelse(k == 1, logTrue, logRec)), hold on, plotPatches2d(rcv)
    plot(m.obs(sta,1)/1e3, 0*sta, 'v', 'Color', st.ink, 'MarkerFaceColor', st.ink, 'MarkerSize', 4)
    title(ifelse(k == 1, 'true', sprintf('recovered (%d KL modes)', r)))
    clim([-1 1]*max(abs(logTrue)))
end
colormap(bluewhitered(256))
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'log_{10}(\eta / 10^{19} Pa s)';
xlabel(tl, 'distance from the trench (km)', 'FontSize', 11, 'Color', st.ink2)
ylabel(tl, 'depth (km)', 'FontSize', 11, 'Color', st.ink2)
st.axes()

st.figure(2, 1000, 330)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
semilogy(0:nit, hist(:,1), 'o-', 'Color', st.series(1,:), 'LineWidth', 2), grid on
yline(1, ':', 'noise level', 'Color', st.ink2)
xlabel('Gauss-Newton iteration'), ylabel('\chi^2 / N'), title('Data misfit')
nexttile, hold on
show = [5 9 13];                                      % stations (of the 21)
for i = 1:numel(show)
    j = show(i);
    plot(tObs/yr, 1e2*obs.d(j,:), '.', 'Color', st.series(i,:), 'HandleVisibility', 'off')
    plot(tObs/yr, 1e2*(res.u(j,:)), '-', 'Color', st.series(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('x = %.0f km', m.obs(sta(j),1)/1e3))
    plot(tObs/yr, 1e2*(J0.u(j,:)), '--', 'Color', st.series(i,:), 'LineWidth', 1, 'HandleVisibility', 'off')
end
grid on, legend('Location', 'best', 'Box', 'off')
xlabel('time after the earthquake (yr)'), ylabel('horizontal displacement (cm)')
title('Fit: start (dashed), final (solid), data (dots)')
st.axes()

function [fwd, balance] = solveAt(xi, balance, eta0, Phi, evl, rcv, shz, locked, sc, T, tObs)
% the forward model at parameters xi: Maxwell mantle, fixed fault
shz.alpha = 1 ./ (eta0 * exp(Phi*xi));
sys = assembleLinearSystem(evl, rcv, shz, locked, sc);
seq = struct('T', T, 'b0', sys.b, 'age', T, 'tEvents', 0, 'bEvents', sys.b, ...
    'tObs', tObs, 'tRef', tObs(1));
fwd = solveLinearSequence(sys, seq, 'Balance', balance);
balance = fwd.balance;
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
