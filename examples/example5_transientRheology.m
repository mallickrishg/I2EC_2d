% EXAMPLE 5: transient mantle rheology (Burgers, Andrade), semi-analytical solver
%
% The cycle of example 2 with three mantle rheologies that share the same
% steady-state viscosity (1e19 Pa s), so the same long-term flow:
%   - Maxwell: the steady-state dashpot only (example 2);
%   - Burgers: plus one Kelvin element (eta_K = 1e18 Pa s, mu_K = 30 GPa,
%     relaxation time eta_K/mu_K of about 1 yr);
%   - Andrade: plus transient creep beta*sigma*t^(1/3), fitted by a
%     Prony series of Kelvin elements (computeAndradePronySeries).
% Kelvin elements are set per cell in shz.etaK and shz.muK; newtonianCycle
% keeps them, and assembleLinearSystem adds their transient strain rates
% to the state, so the periodic cycle is still solved without spin-up.
%
% Run from any folder:  run examples/example5_transientRheology
%
% AUTHORS
%   Rishav Mallick, 2026

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; Vpl = m.Vpl; yr = m.yr;
Trecur = 500*yr;

%% the three rheologies, as Kelvin elements per cell
etaM = 1e19*1e-6;                       % steady-state viscosity (MPa s), as in newtonianCycle
burgers.etaK = 1e18*1e-6 * ones(shz.N, 1);
burgers.muK  = m.mu * ones(shz.N, 1);
% Andrade: transient creep equal to the elastic strain sigma/mu after 1 yr
beta = 1/(m.mu*(1*yr)^(1/3));           % 1/(MPa s^(1/3))
fit = computeAndradePronySeries(beta*ones(shz.N, 1), 0.1*yr, 3*Trecur);
andrade.etaK = fit.etaK;
andrade.muK  = fit.muK;
fprintf('Andrade: %d Kelvin elements, tau = %.2g to %.2g yr, fit error %.1e on [0.1, %g] yr\n', ...
    numel(fit.tau), fit.tau(1)/yr, fit.tau(end)/yr, fit.maxRelErr, 3*Trecur/yr);

rheo = {'Maxwell', struct('etaK', zeros(shz.N, 0), 'muK', zeros(shz.N, 0)); ...
        'Burgers', burgers; 'Andrade', andrade};
tvec = logspace(-2, log10(Trecur/yr), 80)*yr;      % times after the earthquake, up to T
xsite = 200;                                        % surface site (km)
[~, isite] = min(abs(m.obs(:,1)/1e3 - xsite));
% surface velocity from mantle flow alone (the fault's afterslip, much
% faster and larger, is left out to show the mantle's response)
Gm = [m.devl.LO(isite,:,1,1), m.devl.LO(isite,:,1,2)];
w = shz.area/sum(shz.area);

edotMean = zeros(size(rheo, 1), numel(tvec));
vxMantle = zeros(size(rheo, 1), numel(tvec));
for r = 1:size(rheo, 1)
    shz.etaK = rheo{r,2}.etaK; shz.muK = rheo{r,2}.muK;
    tic
    cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);
    sol = cyc.solve(tvec);
    % cycle mean (slip and strain over the whole cycle / T): the Kelvin
    % elements change the timing of the deformation, not its average
    Vm = sol.slip(:,end)/Trecur; em = [sol.e22(:,end); sol.e23(:,end)]/Trecur;
    epl = [shz.e22pl; shz.e23pl];
    fprintf(['%-8s %5d states, cycle-mean departure from the long-term rates: ' ...
        'slip %.1e, strain %.1e  (%.0f s)\n'], rheo{r,1}, size(cyc.A, 1), ...
        max(abs(Vm(~cyc.locked) - Vpl))/Vpl, max(abs(em - epl))/max(abs(epl)), toc);
    % mantle: area-weighted mean strain rate, and its surface velocity
    edot = sqrt(sol.e22dot.^2 + sol.e23dot.^2);
    edot_lt = sqrt(shz.e22pl.^2 + shz.e23pl.^2);       % set by newtonianCycle
    edotMean(r,:) = w' * edot / (w' * edot_lt);
    vxMantle(r,:) = Gm * [sol.e22dot - shz.e22pl; sol.e23dot - shz.e23pl];
end
shz.etaK = zeros(shz.N, 0); shz.muK = zeros(shz.N, 0);   % leave shz Maxwell
fprintf('mean mantle strain rate / long-term at %.2g yr: %s; at %.0f yr: %s\n', ...
    tvec(1)/yr, mat2str(edotMean(:,1)', 3), tvec(end)/yr, mat2str(edotMean(:,end)', 3));

st = exampleStyle();
col = st.series;

%% creep under constant stress: the three rheologies
st.figure(1, 760, 380), hold on
t = logspace(-3, log10(3*Trecur/yr), 300)*yr;
J0 = t/etaM;                                         % steady-state creep
JB = J0 + (1 - exp(-t*burgers.muK(1)/burgers.etaK(1)))/burgers.muK(1);
JA = J0 + beta*t.^(1/3);
JP = J0 + ((1 - exp(-t(:)./fit.tau))*(1./andrade.muK(1,:))')';
plot(t/yr, J0*m.mu, 'Color', col(1,:), 'LineWidth', 2, 'DisplayName', 'Maxwell')
plot(t/yr, JB*m.mu, 'Color', col(2,:), 'LineWidth', 2, 'DisplayName', 'Burgers')
plot(t/yr, JA*m.mu, 'Color', col(3,:), 'LineWidth', 2, 'DisplayName', 'Andrade, \beta t^{1/3}')
plot(t/yr, JP*m.mu, '--', 'Color', st.ink, 'LineWidth', 1.25, ...
    'DisplayName', sprintf('Andrade, Prony series (%d terms)', numel(fit.tau)))
xline([0.1 3*Trecur/yr], ':', 'Color', st.ink2, 'HandleVisibility', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log'), grid on
xlim(t([1 end])/yr), ylim([1e-3 1e3])
legend('Location', 'northwest', 'Box', 'off')
xlabel('time under constant stress (yr)'), ylabel('viscous strain / elastic strain')
title('Creep of one cell (dotted: Prony fit range)')
st.axes()

%% postseismic mantle response
st.figure(2, 1000, 380)
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile, hold on
for r = 1:size(rheo, 1)
    plot(tvec/yr, edotMean(r,:), 'Color', col(r,:), 'LineWidth', 2, 'DisplayName', rheo{r,1})
end
yline(1, ':', 'long-term', 'Color', st.ink2, 'HandleVisibility', 'off', 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom')
set(gca, 'XScale', 'log', 'YScale', 'log'), grid on, xlim(tvec([1 end])/yr)
ylim([0.5 2*max(edotMean(:))])
legend('Location', 'northeast', 'Box', 'off')
ylabel('mean mantle strain rate / long-term')
title('Mantle strain rate (area-weighted mean)')
nexttile, hold on
yline(0, 'Color', st.grid, 'LineWidth', 1, 'HandleVisibility', 'off')
for r = 1:size(rheo, 1)
    plot(tvec/yr, vxMantle(r,:)/Vpl, 'Color', col(r,:), 'LineWidth', 2, 'DisplayName', rheo{r,1})
end
set(gca, 'XScale', 'log'), grid on, xlim(tvec([1 end])/yr)
ylabel('horizontal velocity / plate rate')
title(sprintf('Surface velocity from mantle flow, x = %d km', xsite))
xlabel(tl, 'time after the earthquake (yr)', 'FontSize', 11, 'Color', st.ink2)
st.axes()
