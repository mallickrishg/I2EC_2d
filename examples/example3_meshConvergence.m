% EXAMPLE 3: mesh convergence of the mantle strain rate over a cycle
%
% Solves the periodic cycle of example 2 on three mantle meshes of the
% same geometry and spatial domain (coarse, regular, fine; each halves
% the largest triangle edge of the previous one) and compares snapshots
% of the mantle strain rate at t/T = 0.01, 0.1 and 0.99 after the
% earthquake. The semi-analytical solver is exact in time, so any
% difference between the three solutions comes from the mesh. The fault
% and boundary meshes are the same in all three runs.
%
% Meshes are made with meshing/createShearZoneMesh (PDE Toolbox) and
% their kernels (default: cell-averaged, certified) are cached as HDF5
% files in kernels/meshconvergence/ (git-ignored). The first run computes them,
% about an hour for the fine mesh; later runs take about a minute.
% Run from any folder:  run examples/example3_meshConvergence

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel();
yr = m.yr;

% snapshot times after the earthquake, as a fraction of the cycle
tfrac = [0.01 0.1 0.99];

% mantle meshes: largest triangle edge (km)
res = struct('name', {'coarse', 'regular', 'fine'}, 'Hmax', {40, 20, 10});

cachedir = fullfile(m.root, 'kernels', 'meshconvergence');
if ~exist(cachedir, 'dir'), mkdir(cachedir); end

%% solve the periodic cycle on each mesh
for k = 1:numel(res)
    name = sprintf('shearzone_h%02d', res(k).Hmax);
    if ~isfile(fullfile(cachedir, [name '_vertices.dat']))
        createShearZoneMesh('Hmax', res(k).Hmax, 'OutputFolder', cachedir, 'Name', name);
    end
    mk = exampleModel('ShearZone', fullfile(cachedir, name));   % fresh receivers
    shz = mk.shz;

    kernelfile = fullfile(cachedir, ['evl_' name '.h5']);
    oldfile = fullfile(cachedir, ['evl_' name '.mat']);   % cache of earlier versions
    if isfile(kernelfile)
        evl = loadKernels(kernelfile);
    elseif isfile(oldfile)
        load(oldfile, 'evl');
        saveKernels(kernelfile, evl, 'rcv', mk.rcv, 'shz', shz);
    else
        fprintf('computing kernels for the %s mesh (%d cells)\n', res(k).name, shz.N);
        evl = computeAllStressKernelsBem(mk.rcv, shz, mk.boundary);
        saveKernels(kernelfile, evl, 'rcv', mk.rcv, 'shz', shz);
    end
    assert(size(evl.LL,1) == shz.N, 'cached kernel does not match %s', name);

    cyc = newtonianCycle(mk.rcv, shz, evl, mk.Vpl, yr);
    sol = cyc.solve(tfrac*cyc.Trecur);
    res(k).edot = sqrt(sol.e22dot.^2 + sol.e23dot.^2);
    res(k).shz = shz;
    res(k).rcv = mk.rcv;
    res(k).Trecur = cyc.Trecur;
    fprintf('%-8s Hmax = %2d km: %4d cells\n', res(k).name, res(k).Hmax, shz.N);
end

%% compare: every mesh sampled at the same points
% value of the containing triangle, on a 2 km grid over the domain
[xg, zg] = meshgrid(-100e3:2e3:450e3, -200e3:2e3:-20e3);
for k = 1:numel(res)
    res(k).sample = sampleOnPoints(res(k).shz, res(k).edot, [xg(:) zg(:)]);
end
inside = all(~isnan([res.sample]), 2);   % points inside the mantle on every mesh
fine = res(end).sample(inside,:);

fprintf('\nDifference from the fine mesh, %d points on a 2 km grid\n', sum(inside));
fprintf('%-26s', '');  fprintf('  t/T = %-6.2f', tfrac);  fprintf('\n');
for k = 1:numel(res)-1
    s = res(k).sample(inside,:);
    fprintf('%-8s relative L2       ', res(k).name);
    fprintf('  %9.1f %% ', 100*vecnorm(s - fine)./vecnorm(fine));  fprintf('\n');
    fprintf('%-8s median |log10 ratio|', '');
    fprintf('  %11.3f ', median(abs(log10(s./fine))));  fprintf('\n');
end

%% plot: snapshots, one row per mesh, one column per time
figure(1),clf
tl = tiledlayout(numel(res), numel(tfrac), 'TileSpacing', 'compact');
for k = 1:numel(res)
    for i = 1:numel(tfrac)
        nexttile
        plotshz2d(res(k).shz, res(k).edot(:,i)), hold on, plotpatch2d(res(k).rcv)
        axis tight equal, box on
        clim(10.^[-15 -12])
        set(gca, 'ColorScale', 'log', 'FontSize', 11, 'LineWidth', 1, 'TickDir', 'in')
        if k == 1
            title(sprintf('t/T = %.2f (%.0f yr)', tfrac(i), tfrac(i)*res(k).Trecur/yr))
        end
        if i == 1
            ylabel(sprintf('%s, %d cells\nz (km)', res(k).name, res(k).shz.N))
        end
        if k == numel(res), xlabel('x (km)'), end
    end
end
colormap(turbo(30))
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.Interpreter = 'latex';     % before the string, which is LaTeX
cb.Label.String = 'strain rate $\sqrt{\dot\epsilon_{22}^2 + \dot\epsilon_{23}^2}$ (1/s)';
title(tl, 'Mantle strain rate over the earthquake cycle on three meshes')

%% plot: profiles through all three meshes
profiles = struct('label', {'z = -60 km', 'x = 300 km'}, ...
    'pts', {[linspace(-100e3, 450e3, 551)' -60e3*ones(551,1)], ...
            [300e3*ones(331,1) linspace(-35e3, -200e3, 331)']}, ...
    'horizontal', {true, false});
figure(2),clf
tl = tiledlayout(numel(profiles), numel(tfrac), 'TileSpacing', 'compact');
for p = 1:numel(profiles)
    for i = 1:numel(tfrac)
        nexttile, hold on
        for k = 1:numel(res)
            s = sampleOnPoints(res(k).shz, res(k).edot(:,i), profiles(p).pts);
            if profiles(p).horizontal
                plot(profiles(p).pts(:,1)/1e3, s, 'LineWidth', 1.5)
            else
                plot(s, profiles(p).pts(:,2)/1e3, 'LineWidth', 1.5)
            end
        end
        box on, grid on
        if profiles(p).horizontal
            set(gca, 'YScale', 'log'), xlabel('x (km)'), ylabel('strain rate (1/s)')
        else
            set(gca, 'XScale', 'log'), xlabel('strain rate (1/s)'), ylabel('z (km)')
        end
        set(gca, 'FontSize', 11, 'LineWidth', 1)
        title(sprintf('%s, t/T = %.2f', profiles(p).label, tfrac(i)))
        if p == 1 && i == 1
            legend(arrayfun(@(r) sprintf('%s (%d cells)', r.name, r.shz.N), res, ...
                'UniformOutput', false), 'Location', 'best')
        end
    end
end
title(tl, 'Strain-rate profiles on the three meshes')

function s = sampleOnPoints(shz, values, pts)
% value of the triangle that contains each point (NaN outside the mesh)
TR = triangulation(shz.tri, shz.vert);
id = pointLocation(TR, pts);
s = nan(size(pts,1), size(values,2));
in = ~isnan(id);
s(in,:) = values(id(in),:);
end
