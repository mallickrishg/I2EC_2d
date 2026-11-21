% COUPLEDDISTANCESWEEP  Coupled-operator certification for a fault at distance d from the mantle.
%
% A 100 km fault (dip 14 deg, 2.5 km patches) above a structured mantle
% band (5 km cells) parallel to it at normal distance d. For each d,
% compares the default kernels (KL cell-averaged with 3 points, LK
% collocated at patch centres) with averaged and Betti-consistent
% coupling blocks, and reports the KL/LK mismatch, the largest
% eigenvalue of the symmetric part of the scaled, weighted coupled
% operator, and the growth of K*diag(D) for random positive D.
% See tasks.md, section 0a.
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root, fullfile(root, 'functions'), fullfile(root, 'meshing'));
folder = fullfile(tempdir, 'i2ec_coupled_sweep'); if ~exist(folder, 'dir'), mkdir(folder); end

G = 30e3; nu = 0.25; em = geometry.LDhs(G, nu);
dip = 14; L = 100e3; W0 = 2.5e3; h = 5e3; T = 30e3;
fid = fopen(fullfile(folder, 'fault.seg'), 'w');
fprintf(fid, '# n Vpl x z Width Dip W0 qW\n1 1.0 0.0 1.0 %.6f %.6f %.6f 1 0 0\n', L, dip, W0);
fclose(fid);
fid = fopen(fullfile(folder, 'bnd.seg'), 'w');
fprintf(fid, '# n Vpl x z Width Dip W0 qW bcx bcz\n');
fprintf(fid, '1 0.0 -500e3 0     500e3   90 50e3 1 0 0\n');
fprintf(fid, '2 0.0 -500e3 500e3 1250e3   0 62.5e3 1 0 0\n');
fprintf(fid, '3 0.0  750e3 500e3 500e3  -90 50e3 1 0 0\n');
fclose(fid);
rcv = geometry.receiver(fullfile(folder, 'fault.seg'), em);
bnd = geometry.receiver(fullfile(folder, 'bnd.seg'), em);
N = rcv.N; Wp = rcv.W(:);
sv = [cosd(dip) -sind(dip)];          % down-dip
nv = [-sind(dip) -cosd(dip)];         % into the footwall

ds = [20e3 10e3 5e3 2e3 1e3 0.5e3 0.1e3 0];
rows = {};
for d = ds
    % structured mesh in fault coordinates (s along dip, n normal), one diagonal
    s = 0:h:L; n = d + (0:h:T);
    [S, Nn] = meshgrid(s, n); vert = S(:)*sv + Nn(:)*nv;
    ns = numel(s); nn = numel(n); id = reshape(1:ns*nn, nn, ns);
    tri = [];
    for i = 1:nn-1
        for j = 1:ns-1
            p = id(i,j); q = id(i+1,j); r = id(i,j+1); t = id(i+1,j+1);
            tri = [tri; p q t; p t r]; %#ok<AGROW>
        end
    end
    prefix = fullfile(folder, sprintf('band_%g', d));
    writematrix(vert, [prefix '_vertices.dat']); writematrix(tri, [prefix '_triangulation.dat']);
    shz = geometry.shearZoneReceiver(prefix, em);
    a = shz.area(:); M = shz.N;

    evl3 = computeAllStressKernelsBem(rcv, shz, bnd);                         % default
    evl7 = computeAllStressKernelsBem(rcv, shz, bnd, 'cellaverage', 7);       % finer KL, LL
    LKavg = patchAveragedLK(rcv, shz, bnd, 5);
    LK32 = patchAveragedLK(rcv, shz, bnd, 32);
    KLb = permute(LK32, [2 1 3]) .* (Wp' ./ (2*a));      % Betti: a*KL = (W/2)*LK'
    variants = {'KL 3-pt, LK collocated (default)', evl3; ...
                'KL 3-pt, LK patch-avg 5-pt', setfield(evl3, 'LK', LKavg); ...
                'KL 7-pt, LK patch-avg 5-pt', setfield(evl7, 'LK', LKavg); ...
                'LK patch-avg 32-pt, KL by Betti', setfield(setfield(evl3, 'LK', LK32), 'KL', KLb)}; %#ok<SFLD>
    fprintf('\n=== d = %.1f km (%d cells, %d patches)\n', d/1e3, M, N);
    for v = 1:size(variants, 1)
        r = assess(variants{v,2}, Wp, a);
        fprintf('  %-34s mismatch median %.3f max %.3f | sym max eig %+.2e (%d >= 0) | KD max Re %+.2e (%d > 0)\n', ...
            variants{v,1}, r.medMis, r.maxMis, r.maxEig, r.nPos, r.kdMax, r.kdPos);
        rows(end+1, :) = {d, variants{v,1}, r.medMis, r.maxMis, r.maxEig, r.nPos, r.kdMax, r.kdPos}; %#ok<SAGROW>
    end
end

function r = assess(evl, W, a)
N = numel(W);
K = assembleStressKernel(evl, false(N,1));
w = [W; 2*a; 2*a];
S = w .* K;
dd = 1 ./ sqrt(abs(diag(S)));
St = dd .* S .* dd';
lam = eig((St + St')/2);
C1 = St(1:N, N+1:end); C2 = St(N+1:end, 1:N)';
mis = vecnorm(C1 - C2, 2, 2) ./ vecnorm(C1, 2, 2);
r.medMis = median(mis); r.maxMis = max(mis);
r.maxEig = max(lam); r.nPos = sum(lam >= 0);
% growth of the actual system K*diag(D), D > 0 random (fault and mantle)
rng(2); kd = -inf; np = 0;
for t = 1:10
    D = exp(3*randn(size(K,1), 1));
    mu = eig(K .* D');
    scale = max(abs(mu));
    kd = max(kd, max(real(mu))/scale); np = max(np, sum(real(mu) > 1e-10*scale));
end
r.kdMax = kd; r.kdPos = np;
end

function LK = patchAveragedLK(rcv, shz, bnd, ng)
b = (1:ng-1) ./ sqrt(4*(1:ng-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1)); [s, i] = sort(diag(D)); wq = V(1, i).^2;
bc = computeShzBoundaryCorrection(shz, bnd);
xc0 = rcv.xc; LK = zeros(rcv.N, shz.N, 2);
for q = 1:ng
    rcv.xc = xc0 + (s(q)*rcv.W(:)/2) .* rcv.dv;
    LL = computeShzStressKernelsBem(shz, rcv, bnd, bc);
    LK(:,:,1) = LK(:,:,1) + wq(q)*(LL(:,:,1,1) - LL(:,:,1,3));
    LK(:,:,2) = LK(:,:,2) + wq(q)*LL(:,:,1,2);
end
rcv.xc = xc0;
end
