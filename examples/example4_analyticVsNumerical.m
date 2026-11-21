% EXAMPLE 4: the semi-analytical and the numerical solver on the same cycle
%
% Solves the cycle of example 2 twice:
%   - semi-analytical: periodic state solved directly, evolved exactly
%     (computePeriodicInitialConditionBalanced, analyticSolveIVP);
%   - numerical: the same equations integrated in time by ode45
%     (viscoeq_sequence_simulation, ode_viscoelastic), with the fault law
%     set to linear viscous (rcv.faultLaw = 'linear', done by
%     newtonianCycle) so that both solve the same linear ODE.
% The numerical run starts from the analytic state just before the
% earthquake and integrates one cycle; it must follow the analytic
% solution at every time and end where it started.
%
% Starting the numerical run from the long-term state instead would need
% a spin-up: the slowest mode of this model decays over thousands of
% years (printed below), so it takes many 500-year cycles to become
% periodic. The test TimeIntegrationTest/spunUp* checks the spin-up on a
% smaller model.
% Run from any folder:  run examples/example4_analyticVsNumerical

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel('StressKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl; Vpl = m.Vpl; yr = m.yr;
cyc = newtonianCycle(rcv, shz, evl, Vpl, yr);
locked = cyc.locked; Nc = sum(~locked); M = shz.N;
lam = eig(cyc.A);
fprintf('slowest mode: %.0f yr (rho = %.3f per cycle, %d cycles to spin up to 1e-6)\n', ...
    1/min(abs(real(lam)))/yr, max(abs(exp(lam*cyc.Trecur))), ...
    ceil(log(1e-6)/log(max(abs(exp(lam*cyc.Trecur))))));

%% numerical: one cycle from the analytic state just before the earthquake
ypre = cyc.y0 - cyc.b;                         % [v(creeping); e22dot; e23dot]
Y0 = zeros(rcv.N + 2*M, 1);                    % [v; s22, s23 interleaved]
Y0(~locked) = ypre(1:Nc);
Y0(rcv.N+1:2:end) = ypre(Nc+1:Nc+M) ./ shz.alpha;
Y0(rcv.N+2:2:end) = ypre(Nc+M+1:end) ./ shz.alpha;
tic
[t, Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl), ...
    Y0, rcv, shz, cyc.stress_change, cyc.Trecur);       % earthquake at t = 0
fprintf('numerical cycle: %d time steps, %.1f s\n', numel(t), toc);
[Vn, e22n, e23n] = extract_v_edot_from_statevector(Y, rcv, shz);

%% analytic: the same times
sol = cyc.solve(t');
Va = sol.V'; e22a = sol.e22dot'; e23a = sol.e23dot';

%% compare
relV = max(abs(Vn - Va), [], 'all') / max(abs(Va), [], 'all');
relE = max(abs([e22n e23n] - [e22a e23a]), [], 'all') / max(abs([e22a e23a]), [], 'all');
used = [~locked; true(2*M, 1)];                % locked fault states are not used
closure = max(abs(Y(end,used)' - Y0(used)) ./ max(abs(Y0(used))));
fprintf('largest difference over the cycle: slip rate %.1e, strain rate %.1e (relative)\n', relV, relE);
fprintf('numerical end state vs start state: %.1e (relative)\n', closure);

%% plot: time series at a few fault patches and mantle cells
[~, ip] = min(abs(rcv.xc(:,2) + [32 40 48]*1e3));          % patches at 32, 40, 48 km depth
pts = [150 -60; 250 -60; 300 -150]*1e3;                    % mantle points (m)
ic = arrayfun(@(k) find(inShzCell(shz, pts(k,:)), 1), 1:size(pts,1));
tt = t/cyc.Trecur;
keep = tt > 0;

st = exampleStyle();
st.figure(1, 1000, 420)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile, hold on
for k = 1:numel(ip)
    plot(tt(keep), Va(keep,ip(k))/Vpl, '-', 'Color', st.series(k,:), 'LineWidth', 3, ...
        'DisplayName', sprintf('%.0f km depth', -rcv.xc(ip(k),2)/1e3));
    plot(tt(keep), Vn(keep,ip(k))/Vpl, '--', 'Color', st.ink, 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(nan, nan, '--', 'Color', st.ink, 'LineWidth', 1, 'DisplayName', 'numerical (ode45)')
set(gca, 'XScale', 'log', 'YScale', 'log'), grid on
xlim([1e-6 1]), xlabel('time since the earthquake, t/T'), ylabel('slip rate / plate rate')
title('Fault slip rate')
legend('Location', 'southwest', 'Box', 'off')
nexttile, hold on
ea = sqrt(e22a.^2 + e23a.^2); en = sqrt(e22n.^2 + e23n.^2);
for k = 1:numel(ic)
    plot(tt(keep), ea(keep,ic(k)), '-', 'Color', st.series(k,:), 'LineWidth', 3, ...
        'DisplayName', sprintf('x = %.0f km, %.0f km depth', pts(k,1)/1e3, -pts(k,2)/1e3));
    plot(tt(keep), en(keep,ic(k)), '--', 'Color', st.ink, 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(nan, nan, '--', 'Color', st.ink, 'LineWidth', 1, 'DisplayName', 'numerical (ode45)')
set(gca, 'XScale', 'log', 'YScale', 'log'), grid on
xlim([1e-6 1]), xlabel('time since the earthquake, t/T'), ylabel('strain rate (1/s)')
title('Mantle strain rate')
legend('Location', 'southwest', 'Box', 'off')
st.axes()

function in = inShzCell(shz, p)
% true for the triangle(s) of shz that contain point p
TR = triangulation(shz.tri, shz.vert);
in = false(shz.N, 1);
id = pointLocation(TR, p);
if ~isnan(id), in(id) = true; end
end
