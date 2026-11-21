% COUPLEDZEROGAP  Example model with the fault tip touching the mantle wedge.
%
% In examples/inputs/shearzone the left edge of the mantle wedge is the
% line x = 201 km (z = -35 to -50 km), so the down-dip tip of the
% megathrust (x = 200 km, z = -50 km) is 1 km from the wedge corner. This
% script moves that edge onto the fault tip (gap 0), rebuilds the default
% kernels (3-point, certified), checks the width/area-weighted coupled
% operator, and solves the example 2 Newtonian cycle on it.
% See tasks.md, section 0a.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root, 'examples'));
m = exampleModel();
rcv = m.rcv;
ends = rcv.xc(end,:) + [-1; 1] * (rcv.W(end)/2 * rcv.dv(end,:));
[~, k] = min(ends(:,2)); tip = ends(k,:);                % deeper end of patch 49
fprintf('fault tip at (%.3f, %.3f) km\n', tip/1e3);

V = readmatrix(fullfile(root, 'examples', 'inputs', 'shearzone_vertices.dat'));
T = readmatrix(fullfile(root, 'examples', 'inputs', 'shearzone_triangulation.dat'));
edge = abs(V(:,1) - 201e3) < 1 & V(:,2) >= -50e3 - 1 & V(:,2) <= -35e3 + 1;
fprintf('moving %d wedge-edge vertices from x = 201 km to the tip\n', sum(edge));
V(edge, 1) = tip(1);
[~, corner] = min(abs(V(:,2) - (-50e3)) + 1e9*~edge);
V(corner, :) = tip;

folder = fullfile(tempdir, 'i2ec_zero_gap'); if ~exist(folder, 'dir'), mkdir(folder); end
prefix = fullfile(folder, 'shearzone_zerogap');
writematrix(V, [prefix '_vertices.dat']); writematrix(T, [prefix '_triangulation.dat']);
mk = exampleModel('ShearZone', prefix);
shz = mk.shz; rcv = mk.rcv;
fprintf('smallest cell area %.2f km^2 (all positive: %d)\n', min(shz.area)/1e6, all(shz.area > 0));

evl = computeAllStressKernelsBem(rcv, shz, mk.boundary);      % default: 3-point, certified
fprintf('LL certificate: asymmetry %.3e, max eigenvalue %.3e\n', ...
    evl.certificate.asymmetry, evl.certificate.maxEigenvalue);

N = rcv.N; W = rcv.W(:); a = shz.area(:);
K = assembleStressKernel(evl, false(N,1));
S = [W; 2*a; 2*a] .* K;
dd = 1 ./ sqrt(abs(diag(S)));
St = dd .* S .* dd';
lam = eig((St + St')/2);
C1 = St(1:N, N+1:end); C2 = St(N+1:end, 1:N)';
mis = vecnorm(C1 - C2, 2, 2) ./ vecnorm(C1, 2, 2);
fprintf('coupled operator (scaled): sym part max eig %.3e, %d >= 0\n', max(lam), sum(lam >= 0));
fprintf('KL/LK mismatch per patch: median %.3f; patches %d-%d: %s\n', median(mis), N-3, N, ...
    mat2str(round(mis(end-3:end)', 3)));

rng(1); worst = -inf;
for trial = 1:20
    D = exp(4*randn(size(K,1), 1));
    mu = eig(K .* D');
    worst = max(worst, max(real(mu))/max(abs(mu)));
end
fprintf('K*diag(D), 20 random D: max Re/max|.| = %.3e\n', worst);

cyc = newtonianCycle(rcv, shz, evl, mk.Vpl, mk.yr);
sol = cyc.solve([0.01 0.5 0.99]*cyc.Trecur);
fprintf('example 2 cycle solved: max |e22dot| %.3e 1/s, all finite: %d\n', ...
    max(abs(sol.e22dot(:))), all(isfinite(sol.e22dot(:))));
