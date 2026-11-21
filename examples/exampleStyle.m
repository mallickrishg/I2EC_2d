function s = exampleStyle()
% EXAMPLESTYLE  Colours and axes style shared by the example figures.
%
%   s = exampleStyle()
%
% OUTPUTS
%   s   struct with fields
%       .series   [3 x 3] colours for unordered series (blue, orange, aqua)
%       .times    [3 x 3] colours for three ordered times, dark (early) to
%                 light (late)
%       .blue     [256 x 3] sequential colormap, light to dark blue
%       .orange   [256 x 3] sequential colormap, light to dark orange
%       .ink, .ink2, .grid   text, axis and gridline colours
%       .axes     @(ax) applies the shared axes style to ax (or to every
%                 axes of the current figure when called without input)
%       .map      @(ax) the style for a cross-section map: equal axes,
%                 faint cell edges, no grid, depth labelled positive down
%       .figure   @(n, w, h) opens figure n, cleared, w x h pixels
%
% NOTES
%   Colours follow a validated data-visualisation palette: one hue, light
%   to dark, for magnitudes; fixed categorical slots for unordered series.

hex = @(h) sscanf(h(2:end), '%2x%2x%2x', [1 3])/255;
ramp = @(hs) interp1(linspace(0, 1, numel(hs)), cell2mat(cellfun(hex, hs', 'UniformOutput', false)), ...
    linspace(0, 1, 256));

s.series = [hex('#2a78d6'); hex('#eb6834'); hex('#1baf7a')];
s.times  = [hex('#0d366b'); hex('#2a78d6'); hex('#86b6ef')];
s.blue   = ramp({'#f2f7fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'});
s.orange = ramp({'#fdf1ea', '#fbd5c0', '#f6ae8a', '#f08655', '#eb6834', '#c5501f', '#933a14', '#5e240b'});
s.ink  = hex('#0b0b0b');
s.ink2 = hex('#52514e');
s.grid = hex('#e1e0d9');
s.axes = @(varargin) styleAxes(s, varargin{:});
s.map  = @(varargin) styleMap(s, varargin{:});
s.figure = @(n, w, h) openFigure(n, w, h);
end

function styleAxes(s, ax)
if nargin < 2, ax = findobj(gcf, 'Type', 'axes'); end
set(ax, 'FontName', 'Helvetica', 'FontSize', 11, 'LineWidth', 0.75, ...
    'TickDir', 'out', 'TickLength', [0.008 0.008], 'Box', 'off', ...
    'XColor', s.ink2, 'YColor', s.ink2, 'GridColor', s.grid, 'GridAlpha', 1, ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'TitleFontWeight', 'normal', 'TitleFontSizeMultiplier', 1.1, 'Layer', 'top');
for a = ax(:)'
    a.Title.Color = s.ink;
end
cb = findobj(gcf, 'Type', 'colorbar');
set(cb, 'FontName', 'Helvetica', 'FontSize', 10, 'Color', s.ink2, 'TickDirection', 'out', ...
    'LineWidth', 0.75, 'Box', 'off');
end

function styleMap(s, ax)
if nargin < 2, ax = findobj(gcf, 'Type', 'axes'); end
styleAxes(s, ax)
for a = ax(:)'
    axis(a, 'equal', 'tight')
    set(findobj(a, 'Type', 'patch'), 'EdgeColor', 'k', 'EdgeAlpha', 0.07, 'LineWidth', 0.25)
    set(findobj(a, 'Type', 'line'), 'Color', s.ink, 'LineWidth', 2)
    a.Box = 'on'; a.XGrid = 'off'; a.YGrid = 'off';
    a.YTickLabel = compose('%g', abs(a.YTick));      % depth, positive down
end
end

function openFigure(n, w, h)
f = figure(n); clf(f)
f.Units = 'pixels';
f.Position(3:4) = [w h];
end
