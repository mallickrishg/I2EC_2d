% EXAMPLE 1: earthquake cycles with a power-law (n = 3) mantle, numerical solver
%
% The geometry, locking and earthquake of example 2, with a power-law
% mantle, edot = alpha * tau^(n-1) * s with n = 3 everywhere, and a
% rate-strengthening creeping fault, d log(v)/dt = (dtau/dt)/(a-b)sigma.
% The power-law coefficient alpha is uniform, so the mantle is one
% power-law material and the n = 3 corner flow (getStrainratesLongterm)
% is its long-term (plate-rate) solution. alpha gives an effective
% viscosity of 1e19 Pa s, the viscosity of example 2, at the median
% long-term strain rate.
%
% The system is nonlinear, so it is integrated in time (ode45 in
% runImposedEarthquakeCycles) from the long-term state through several
% cycles until it is close to periodic; the last cycle is plotted.
% Run from any folder:  run examples/example1_numericalPowerLaw
%
% AUTHORS
%   Sharadha Sathiakumar & Rishav Mallick, EOS 2023

clear
addpath(fileparts(mfilename('fullpath')))
% meshes and kernels: the default kernels (3-point cell average,
% certified free of growing modes); see exampleModel for the earlier
% centroid kernels and their options
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl; Vpl = m.Vpl; yr = m.yr;

Trecur = 500*yr;      % earthquake cycle (s)
Ncycles = 15;         % cycles to integrate (spin-up + the plotted one)
tau_max = 5;          % largest coseismic stress step on the creeping fault (MPa)

%% rheology
% fault: rate-strengthening, (a-b)sigma = 0.2 MPa; locked above 30 km
rcv.faultLaw = 'ratestrengthening';
rcv.Asigma = 0.2 .* ones(rcv.N, 1);
locked = abs(rcv.xc(:,2)) < 30e3;
rcv.pinnedPosition = locked;
rcv.Vpl = Vpl .* ones(rcv.N, 1);

% long-term strain rates (corner flow)
[e22_lt, e23_lt] = getStrainratesLongterm(shz, rcv.dip(1)*pi/180, [0, 20e3], [-140e3, 35e3], 3);
shz.e22pl = e22_lt .* Vpl;
shz.e23pl = -e23_lt .* Vpl;
edot_lt = sqrt(shz.e22pl.^2 + shz.e23pl.^2);

% mantle: n = 3 with one alpha everywhere (a uniform power-law material,
% for which the n = 3 corner flow above is the long-term solution). alpha
% gives an effective viscosity tau/edot of 1e19 Pa s at the median
% long-term strain rate; cells deforming faster in the long term are
% weaker, slower ones stronger (eta_eff ~ edot^(1/n - 1)).
shz.n = 3 .* ones(shz.N, 1);
eta_ref = 1e19*1e-6;                             % MPa s
edot_ref = median(edot_lt);                      % 1/s
shz.alpha = edot_ref ./ (eta_ref*edot_ref).^3 .* ones(shz.N, 1);   % 1/(MPa^3 s)
eta_eff = edot_lt.^(1/3 - 1) .* shz.alpha(1).^(-1/3) * 1e6;          % Pa s
fprintf('long-term effective viscosity: %.2g to %.2g Pa s (1e19 at %.2g 1/s)\n', ...
    min(eta_eff), max(eta_eff), edot_ref);

%% earthquake: the slip deficit of one cycle on the locked zone
slip = Trecur*Vpl .* locked;
dtau = min(evl.KK*slip, tau_max);                % cap the step next to the locked zone
dtau(locked) = 0;
stress_change = struct('Nevents', 1, 'Timing', 1, 'dtau', dtau, ...   % at t = 1 s
    'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);

%% integrate
tic
[t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, stress_change, Ncycles, Trecur);
toc
edot = sqrt(e22dot.^2 + e23dot.^2);
% a periodic cycle ends where it started (both rows are before the
% earthquake): relative change of each cell's strain rate and each
% creeping patch's slip rate over the last cycle
dE = abs(edot(end,:) - edot(1,:)) ./ edot(1,:);
dV = abs(V(end,~locked) - V(1,~locked)) ./ V(1,~locked);
fprintf('change over the last cycle: strain rate median %.1e, max %.1e; slip rate max %.1e (relative, per cell/patch)\n', ...
    median(dE), max(dE), max(dV));

%% fault slip rate through the cycle
st = exampleStyle();
Vplot = V'/Vpl;  Vplot(locked,:) = NaN;          % locked zone: blank
st.figure(1, 760, 420)
pcolor(t/Trecur, rcv.xc(:,1)/1e3, Vplot), shading interp, hold on
text(1e-2, mean(rcv.xc(locked,1))/1e3, 'locked zone (no slip)', 'Color', st.ink2, ...
    'HorizontalAlignment', 'center', 'FontSize', 11)
xlabel('time since the earthquake, t/T'), ylabel('distance from the trench (km)')
cb = colorbar; cb.Label.String = 'slip rate / plate rate';
clim(10.^[-1 2]), colormap(st.orange)
set(gca, 'ColorScale', 'log', 'XScale', 'log', 'YDir', 'reverse')
xlim([1e-4 1]), title('Slip rate on the megathrust')
st.axes()

%% mantle strain rate snapshots
tfrac = [0.01 0.1 0.99];
idx = arrayfun(@(f) find(t >= f*Trecur, 1), tfrac);
st.figure(2, 1000, 470)
tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plotshz2d(shz, edot_lt), hold on, plotpatch2d(rcv)
title('long-term')
for i = 1:numel(tfrac)
    nexttile
    plotshz2d(shz, edot(idx(i),:)'), hold on, plotpatch2d(rcv)
    title(sprintf('%.0f yr after the earthquake (t/T = %.2f)', t(idx(i))/yr, t(idx(i))/Trecur))
end
ax = findobj(gcf, 'Type', 'axes');
set(ax, 'ColorScale', 'log', 'CLim', 10.^[-15 -12])
xlabel(tl, 'distance from the trench (km)', 'FontSize', 11, 'Color', st.ink2)
ylabel(tl, 'depth (km)', 'FontSize', 11, 'Color', st.ink2)
colormap(st.blue)
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'strain rate (1/s)';
st.map()

%% surface velocity
G = @(j) [m.devl.KO(:,:,j), m.devl.LO(:,:,j,1), m.devl.LO(:,:,j,2)];
dstate = [V(idx,:)'; e22dot(idx,:)'; e23dot(idx,:)'] - [rcv.Vpl; shz.e22pl; shz.e23pl];
vx = G(1)*dstate - m.Ghx*(m.hinge.Vpl*Vpl);
vz = G(2)*dstate - m.Ghz*(m.hinge.Vpl*Vpl);
xo = m.obs(:,1)/1e3;
st.figure(3, 760, 520)
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile, hold on
yline(0, 'Color', st.grid, 'LineWidth', 1, 'HandleVisibility', 'off')
for i = 1:numel(idx)
    plot(xo, vx(:,i)/Vpl + (xo <= 0), 'Color', st.times(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%.0f yr (t/T = %.2f)', t(idx(i))/yr, t(idx(i))/Trecur))
end
grid on, ylabel('horizontal / plate rate')
legend('Location', 'northeast', 'Box', 'off')
title('Surface velocity after the earthquake')
nexttile, hold on
yline(0, 'Color', st.grid, 'LineWidth', 1)
for i = 1:numel(idx)
    plot(xo, vz(:,i)/Vpl, 'Color', st.times(i,:), 'LineWidth', 2)
end
grid on, ylabel('vertical / plate rate'), xlabel('distance from the trench (km)')
st.axes()
