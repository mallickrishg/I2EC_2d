function [V, Vx, Vy, lambda] = computeKLmodesMesh(mesh, L, m)
% COMPUTE_EIGMODES_MESH  Build kernel eigenbasis and derivative operators.
%   [V, Vx, Vy, lambda] = compute_eigmodes_mesh(mesh, L)
%   [V, Vx, Vy, lambda] = compute_eigmodes_mesh(mesh, L, m)
%
% INPUT
%   mesh.xc          : (N x 2) centroid coords ([x,z])
%   L                : kernel length scale
%   m (optional)     : number of modes to return (default: all)
%
% OUTPUT
%   V   : (N x m) mode matrix (columns orthonormal)
%   Vx  : (N x m) ∂/∂x of modes
%   Vy  : (N x m) ∂/∂y of modes
%   lambda : (m x 1) singular values ( = eigenvalues for symmetric PD K )
%
% NOTES:
%   - Kernel: k(r) = exp(-(r/L)^2)
%   - If you want full set of N modes omit m or set m = N (but consider memory/time).

if nargin < 3 || isempty(m)
    m = []; % empty -> compute full
end

xc = mesh.xc(:,1);
yc = mesh.xc(:,2);
N  = numel(xc);

% build pairwise differences and squared distances
dx = xc - xc.';      % (N x N)
dy = yc - yc.';
r2 = dx.^2 + dy.^2;

% kernel matrix
K = exp( - r2 / L^2 );

% choose decomposition method
if isempty(m) || m >= N
    % full decomposition: use svd for numerical robustness
    % Note: for symmetric PD K, U*S*V' -> U==V and S diag = eigenvalues
    [U,S] = svd(K, 'econ');   % econ if tall wide; returns U (N x r), S (r x r)
    svals = diag(S);
    V = U; lambda = svals;
else
    % truncated: compute top-m singular vectors/values
    opts.disp = 0;
    % svds returns U,S,V. For symmetric K, U==V. use svds for dense large K.
    [U,S] = svds(K, m, 'largest', opts);
    svals = diag(S);
    V = U; lambda = svals;
end

% precompute derivative kernel matrices (N x N)
r = sqrt(r2);
epsr = 1e-14;
nz = r > epsr;

% radial derivative
k = K;
k_r = zeros(N,N);
k_r(nz) = -(2/L^2) .* r(nz) .* k(nz);

dkdx = zeros(N,N);
dkdy = zeros(N,N);
dkdx(nz) = k_r(nz) .* (dx(nz) ./ r(nz));
dkdy(nz) = k_r(nz) .* (dy(nz) ./ r(nz));

% allocate outputs (N x m)
mout = size(V,2);
Vx = zeros(N, mout);
Vy = zeros(N, mout);

% compute derivative columns: use vector product dkdx * v_k, scaled by 1/lambda_k
for kk = 1:mout
    vk = V(:,kk);
    lam = lambda(kk);
    % avoid division by zero
    if lam == 0
        Vx(:,kk) = 0;
        Vy(:,kk) = 0;
    else
        Vx(:,kk) = (1/lam) * (dkdx * vk);
        Vy(:,kk) = (1/lam) * (dkdy * vk);
    end
end

% (optional) sort modes by descending lambda (largest first)
[lambda, ii] = sort(lambda,'descend');
V  = V(:,ii);
Vx = Vx(:,ii);
Vy = Vy(:,ii);

end
