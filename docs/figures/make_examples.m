% MAKE_EXAMPLES  Figures and console output of the examples for the docs.
%
% Runs the scripts in examples/ and writes every figure they draw to
% docs/_static/examples/<script>_<n>.png (white background, 120 dpi) and
% their console output to docs/_build/examples/<script>.txt.
% example3_meshConvergence is run only when its kernel cache in
% examples/kernels/meshconvergence/ is complete (otherwise its first run builds them).
%
% Run from any folder:
%   run docs/figures/make_examples.m

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
outdir = fullfile(root, 'docs', '_static', 'examples');
logdir = fullfile(root, 'docs', '_build', 'examples');
if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(logdir, 'dir'), mkdir(logdir); end

scripts = {'plotLongtermVelocity', 'example2_analyticNewtonian', ...
           'example4_analyticVsNumerical', 'example1_numericalPowerLaw'};
cache = fullfile(root, 'examples', 'kernels', 'meshconvergence');
if all(cellfun(@(h) exist(fullfile(cache, sprintf('evl_shearzone_h%d.h5', h)), 'file') > 0, {10, 20, 40}))
    scripts{end+1} = 'example3_meshConvergence';
end

set(groot, 'defaultFigureVisible', 'off', 'defaultFigureColor', 'w');
for s = 1:numel(scripts)
    name = scripts{s};
    close all
    fprintf('running %s\n', name);
    logfile = fullfile(logdir, [name '.txt']);
    if exist(logfile, 'file'), delete(logfile); end
    diary(logfile)
    runExample(fullfile(root, 'examples', [name '.m']));
    diary off
    figs = findobj(groot, 'Type', 'figure');
    [~, order] = sort([figs.Number]);
    figs = figs(order);
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(outdir, sprintf('%s_%d.png', name, figs(k).Number)), ...
            'Resolution', 120, 'BackgroundColor', 'white');
    end
end
set(groot, 'defaultFigureVisible', 'remove', 'defaultFigureColor', 'remove');
close all

function runExample(file)
% run in its own workspace: the examples start with "clear"
run(file);
end
