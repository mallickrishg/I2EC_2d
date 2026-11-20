% MAKE_OVERVIEW  Overview figure of the example I2EC model geometry.
%
% Draws the triangulated mantle (oceanic mantle and mantle wedge) and the
% megathrust with its locked zone from the example inputs, and writes
% docs/_static/overview.png (white background, 150 dpi). The model boundary
% (inputs/boundary2d.seg) lies far outside this window and is not drawn.
%
% Run from any folder:
%   run docs/figures/make_overview.m

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(root, fullfile(root, 'functions'));
outfile = fullfile(root, 'docs', '_static', 'overview.png');

% ---- geometry (as in examples/exampleModel.m and newtonianCycle.m)
earthModel = geometry.LDhs(30e3, 0.25);
rcv      = geometry.receiver(fullfile(root, 'inputs', 'megathrust2d.seg'), earthModel);
shz      = geometry.shearZoneReceiver(fullfile(root, 'inputs', 'shearzone'), earthModel);

locked  = abs(rcv.xc(:,2)) < 30e3;                          % locked zone: 0-30 km depth
oceanic = shz.xc(:,1) < -shz.xc(:,2) / tand(rcv.dip(1));    % below the slab

% ---- colours
c_ocean = [0.62 0.78 0.92];
c_wedge = [0.98 0.80 0.60];
c_fault = [0.15 0.15 0.15];
c_lock  = [0.80 0.10 0.10];

km = 1e-3;
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 9], 'Visible', 'off');
ax = axes(fig, 'Position', [0.08 0.14 0.90 0.80]);
hold(ax, 'on');

% triangulated mantle
X = [shz.A(:,1), shz.B(:,1), shz.C(:,1)]' * km;
Z = [shz.A(:,2), shz.B(:,2), shz.C(:,2)]' * km;
patch(ax, X(:, oceanic),  Z(:, oceanic),  c_ocean, 'EdgeColor', [1 1 1]*0.55, 'LineWidth', 0.25);
patch(ax, X(:, ~oceanic), Z(:, ~oceanic), c_wedge, 'EdgeColor', [1 1 1]*0.55, 'LineWidth', 0.25);

% megathrust: creeping part, then locked part on top
[xs, zs] = segments(rcv, ~locked);
h_creep = plot(ax, xs * km, zs * km, '-', 'Color', c_fault, 'LineWidth', 2);
[xs, zs] = segments(rcv, locked);
h_lock = plot(ax, xs * km, zs * km, '-', 'Color', c_lock, 'LineWidth', 3.5);

% free surface
plot(ax, [-150 550], [0 0], 'k-', 'LineWidth', 1);

% labels, placed in empty parts of the panel
slope = atand(tand(rcv.dip(1)));                   % 1:1 aspect: on-screen angle = dip
text(ax, 330, -40, {'mantle wedge', '(continental)'}, 'FontSize', 9, 'HorizontalAlignment', 'center');
text(ax, 100, -150, 'oceanic mantle', 'FontSize', 9, 'HorizontalAlignment', 'center');
text(ax, 330, -92, 'subducting slab', 'FontSize', 9, 'Color', [0.3 0.3 0.3], ...
    'Rotation', -slope, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(ax, 545, -4, 'free surface', 'FontSize', 8, 'Color', [0.3 0.3 0.3], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
legend(ax, [h_lock, h_creep], {'locked megathrust (0–30 km)', 'creeping megathrust'}, ...
    'Location', 'south', 'Orientation', 'horizontal', 'FontSize', 8, 'Box', 'off');

xlim(ax, [-150 550]);
ylim(ax, [-265 5]);
daspect(ax, [1 1 1]);
box(ax, 'on');
set(ax, 'FontSize', 9, 'LineWidth', 0.8, 'TickDir', 'out', 'Layer', 'top');
xlabel(ax, 'Distance from trench (km)');
ylabel(ax, 'Depth (km)');

exportgraphics(fig, outfile, 'Resolution', 150, 'BackgroundColor', 'white');
fprintf('wrote %s\n', outfile);
close(fig);

% ---------------------------------------------------------------------------
function [xs, zs] = segments(p, sel)
% End points of the selected patches as NaN-separated polylines, (x, z) in m.
% p.x(:,2) is already negative below the surface.
k = find(sel(:))';
xs = [p.x(k,1)'; p.x(k,1)' + p.W(k)' .* cosd(p.dip(k)'); nan(1, numel(k))];
zs = [p.x(k,2)'; p.x(k,2)' - p.W(k)' .* sind(p.dip(k)'); nan(1, numel(k))];
xs = xs(:); zs = zs(:);
end
