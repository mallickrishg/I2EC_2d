% MAKE_OVERVIEW  Overview figure of the example I2EC model geometry.
%
% Draws the triangulated mantle (oceanic mantle and mantle wedge), the
% megathrust with its locked zone and the kinematic hinge fault from the
% example inputs. Writes
% docs/_static/overview.png (white background, 150 dpi). The model boundary
% (examples/inputs/boundary2d.seg) lies far outside this window and is not
% drawn.
%
% Run from any folder:
%   run docs/figures/make_overview.m

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(fullfile(root, 'examples'));
outfile = fullfile(root, 'docs', '_static', 'overview.png');

m = exampleModel();                          % geometry as in the examples
rcv = m.rcv; shz = m.shz; hinge = m.hinge;
st = exampleStyle();
locked  = abs(rcv.xc(:,2)) < 30e3;                          % as in the examples
oceanic = shz.xc(:,1) < -shz.xc(:,2) / tand(rcv.dip(1));    % below the slab

c_ocean = st.blue(64,:);
c_wedge = st.orange(56,:);
c_lock  = st.series(2,:);

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [50 50 1100 500], 'Visible', 'off');
ax = axes(fig, 'Position', [0.07 0.13 0.91 0.84]);
drawModel(ax, shz, oceanic, rcv, locked, hinge, c_ocean, c_wedge, c_lock, st, 0.25);
h = ax.UserData;
plot(ax, [-150 550], [0 0], '-', 'Color', st.ink2, 'LineWidth', 0.75);   % free surface

% labels
slope = rcv.dip(1);                          % 1:1 aspect: on-screen angle = dip
text(ax, 330, -55, {'mantle wedge', '(continental)'}, 'FontSize', 10, 'Color', st.ink, ...
    'HorizontalAlignment', 'center');
text(ax, 120, -160, 'oceanic mantle', 'FontSize', 10, 'Color', st.ink, 'HorizontalAlignment', 'center');
text(ax, 320, -97, 'subducting slab (elastic)', 'FontSize', 9, 'Color', st.ink2, ...
    'Rotation', -slope, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(ax, 545, -4, 'free surface', 'FontSize', 9, 'Color', st.ink2, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax, -4, -10, 'hinge', 'FontSize', 9, 'Color', st.ink2, 'HorizontalAlignment', 'right');

xlim(ax, [-150 550]); ylim(ax, [-262 5]); daspect(ax, [1 1 1]);
st.axes(ax); ax.Box = 'off'; ax.XGrid = 'off'; ax.YGrid = 'off';
ax.YTickLabel = compose('%g', abs(ax.YTick));
xlabel(ax, 'distance from the trench (km)'); ylabel(ax, 'depth (km)');
legend(ax, h, {'locked megathrust (0–30 km)', 'creeping megathrust', 'hinge fault'}, ...
    'Location', 'southwest', 'FontSize', 9, 'Box', 'off', 'TextColor', st.ink);

exportgraphics(fig, outfile, 'Resolution', 150, 'BackgroundColor', 'white');
fprintf('wrote %s\n', outfile);
close(fig);

% ---------------------------------------------------------------------------
function drawModel(ax, shz, oceanic, rcv, locked, hinge, c_ocean, c_wedge, c_lock, st, lw)
hold(ax, 'on');
X = [shz.A(:,1), shz.B(:,1), shz.C(:,1)]' / 1e3;
Z = [shz.A(:,2), shz.B(:,2), shz.C(:,2)]' / 1e3;
patch(ax, X(:, oceanic),  Z(:, oceanic),  c_ocean, 'EdgeColor', 'w', 'LineWidth', lw);
patch(ax, X(:, ~oceanic), Z(:, ~oceanic), c_wedge, 'EdgeColor', 'w', 'LineWidth', lw);
[xs, zs] = segments(rcv, ~locked);
h2 = plot(ax, xs/1e3, zs/1e3, '-', 'Color', st.ink, 'LineWidth', 2);
[xs, zs] = segments(rcv, locked);
h1 = plot(ax, xs/1e3, zs/1e3, '-', 'Color', c_lock, 'LineWidth', 3.5);
[xs, zs] = segments(hinge, true(hinge.N, 1));
h3 = plot(ax, xs/1e3, zs/1e3, ':', 'Color', st.ink, 'LineWidth', 1.5);
set([h1 h2 h3], 'Tag', 'legend');
ax.UserData = [h1 h2 h3];
end

function [xs, zs] = segments(p, sel)
% End points of the selected patches as NaN-separated polylines, (x, z) in m.
% p.x(:,2) is already negative below the surface.
k = find(sel(:))';
xs = [p.x(k,1)'; p.x(k,1)' + p.W(k)' .* cosd(p.dip(k)'); nan(1, numel(k))];
zs = [p.x(k,2)'; p.x(k,2)' - p.W(k)' .* sind(p.dip(k)'); nan(1, numel(k))];
xs = xs(:); zs = zs(:);
end
