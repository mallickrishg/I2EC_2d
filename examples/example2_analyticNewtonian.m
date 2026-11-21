% EXAMPLE 2: periodic earthquake cycle, Newtonian mantle, semi-analytical solver
%
% A subduction zone with a megathrust locked above 30 km depth, a linear
% viscous creeping fault below it, and a Newtonian mantle (n = 1) of
% uniform viscosity 1e19 Pa s. An earthquake every 500 years releases the
% slip deficit of the locked zone. The cycle is linear, so its periodic
% state is solved for directly (no spin-up) and evaluated exactly at any
% time (computePeriodicInitialConditionBalanced, analyticSolveIVP).
%
% The model is set up by newtonianCycle, which examples 3 and 4 share.
% Run from any folder:  run examples/example2_analyticNewtonian
%
% AUTHORS
%   Rishav Mallick, JPL 2024; Sharadha Sathiakumar, EOS 2026

clear
addpath(fileparts(mfilename('fullpath')))
% meshes and kernels: the default kernels (3-point cell average,
% certified free of growing modes); see exampleModel for the earlier
% centroid kernels and their options
m = exampleModel('StressKernels', true, 'DisplacementKernels', true);
rcv = m.rcv; shz = m.shz; Vpl = m.Vpl; yr = m.yr;

cyc = newtonianCycle(rcv, shz, m.evl, Vpl, yr);
fprintf('slowest mode decays over %.0f yr, fastest over %.2g yr (cycle: %.0f yr)\n', ...
    1/min(abs(real(eig(cyc.A))))/yr, 1/max(abs(real(eig(cyc.A))))/yr, cyc.Trecur/yr);

% snapshots after the earthquake, as a fraction of the cycle
tfrac = [0.01 0.1 0.99];
sol = cyc.solve(tfrac*cyc.Trecur);
edot = sqrt(sol.e22dot.^2 + sol.e23dot.^2);
edot_lt = sqrt(shz.e22pl.^2 + shz.e23pl.^2);

%% mantle strain rate
st = exampleStyle();
st.figure(1, 1000, 470)
tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plotshz2d(shz, edot_lt), hold on, plotpatch2d(rcv)
title('long-term')
for i = 1:numel(tfrac)
    nexttile
    plotshz2d(shz, edot(:,i)), hold on, plotpatch2d(rcv)
    title(sprintf('%.0f yr after the earthquake (t/T = %.2f)', tfrac(i)*cyc.Trecur/yr, tfrac(i)))
end
set(findobj(gcf, 'Type', 'axes'), 'ColorScale', 'log', 'CLim', 10.^[-15 -12])
xlabel(tl, 'distance from the trench (km)', 'FontSize', 11, 'Color', st.ink2)
ylabel(tl, 'depth (km)', 'FontSize', 11, 'Color', st.ink2)
colormap(st.blue)
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'strain rate (1/s)';
st.map()

%% fault slip rate (creeping part; the locked zone has v = 0)
creep = ~cyc.locked;
st.figure(2, 760, 380), hold on
for i = 1:numel(tfrac)
    plot(rcv.xc(creep,1)/1e3, sol.V(creep,i)/Vpl, 'Color', st.times(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%.0f yr (t/T = %.2f)', tfrac(i)*cyc.Trecur/yr, tfrac(i)))
end
yline(1, ':', 'plate rate', 'Color', st.ink2, 'HandleVisibility', 'off', 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom')
grid on
legend('Location', 'southeast', 'Box', 'off')
xlabel('distance from the trench (km)'), ylabel('slip rate / plate rate')
title('Slip rate on the creeping megathrust')
set(gca, 'YScale', 'log'), ylim padded
st.axes()

%% surface velocity
% perturbation from the long-term state, plus the kinematic hinge
G = @(j) [m.devl.KO(:,:,j), m.devl.LO(:,:,j,1), m.devl.LO(:,:,j,2)];
dstate = [sol.V; sol.e22dot; sol.e23dot] - [rcv.Vpl; shz.e22pl; shz.e23pl];
vx = G(1)*dstate - m.Ghx*(m.hinge.Vpl*Vpl);
vz = G(2)*dstate - m.Ghz*(m.hinge.Vpl*Vpl);
xo = m.obs(:,1)/1e3;

st.figure(3, 760, 520)
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile, hold on
yline(0, 'Color', st.grid, 'LineWidth', 1, 'HandleVisibility', 'off')
for i = 1:numel(tfrac)
    plot(xo, vx(:,i)/Vpl + (xo <= 0), 'Color', st.times(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%.0f yr (t/T = %.2f)', tfrac(i)*cyc.Trecur/yr, tfrac(i)))
end
grid on, ylabel('horizontal / plate rate')
legend('Location', 'northeast', 'Box', 'off')
title('Surface velocity after the earthquake')
nexttile, hold on
yline(0, 'Color', st.grid, 'LineWidth', 1)
for i = 1:numel(tfrac)
    plot(xo, vz(:,i)/Vpl, 'Color', st.times(i,:), 'LineWidth', 2)
end
grid on, ylabel('vertical / plate rate'), xlabel('distance from the trench (km)')
st.axes()
