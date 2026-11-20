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
figure(1),clf
tiledlayout(2, 2, 'TileSpacing', 'compact')
nexttile
plotshz2d(shz, edot_lt), hold on, plotpatch2d(rcv)
title('long-term ($t \to \infty$)', 'Interpreter', 'latex')
for i = 1:numel(tfrac)
    nexttile
    plotshz2d(shz, edot(:,i)), hold on, plotpatch2d(rcv)
    title(sprintf('t/T = %.2f (%.0f yr after the earthquake)', tfrac(i), tfrac(i)*cyc.Trecur/yr))
end
for ax = findobj(gcf, 'Type', 'axes')'
    axis(ax, 'tight', 'equal'), box(ax, 'on'), clim(ax, 10.^[-15 -12])
    set(ax, 'ColorScale', 'log', 'FontSize', 12, 'LineWidth', 1)
    xlabel(ax, 'x (km)'), ylabel(ax, 'z (km)')
end
colormap(turbo(30))
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'strain rate (1/s)';

%% fault slip rate
figure(2),clf
plot(rcv.xc(:,1)/1e3, sol.V/Vpl, 'LineWidth', 2)
grid on, box on
legend(compose('t/T = %.2f', tfrac), 'Location', 'best')
xlabel('x (km)'), ylabel('v/v_{pl}')
title('Slip rate on the megathrust (0 on the locked zone)')
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1)

%% surface velocity
% perturbation from the long-term state, plus the kinematic hinge
G = @(j) [m.devl.KO(:,:,j), m.devl.LO(:,:,j,1), m.devl.LO(:,:,j,2)];
dstate = [sol.V; sol.e22dot; sol.e23dot] - [rcv.Vpl; shz.e22pl; shz.e23pl];
vx = G(1)*dstate - m.Ghx*(m.hinge.Vpl*Vpl);
vz = G(2)*dstate - m.Ghz*(m.hinge.Vpl*Vpl);

figure(3),clf
tiledlayout(2, 1, 'TileSpacing', 'compact')
nexttile
plot(m.obs(:,1)/1e3, vx/Vpl + (m.obs(:,1) <= 0), 'LineWidth', 2)
grid on, box on, ylabel('v_x/v_{pl}')
legend(compose('t/T = %.2f', tfrac), 'Location', 'best')
title('Surface velocity')
nexttile
plot(m.obs(:,1)/1e3, vz/Vpl, 'LineWidth', 2)
grid on, box on, ylabel('v_z/v_{pl}'), xlabel('distance from the trench (km)')
set(findobj(gcf, 'Type', 'axes'), 'FontSize', 12, 'LineWidth', 1)
