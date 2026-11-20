function [V, Vx, Vy, lambda] = computeKLmodesMesh(mesh, L, m)
% COMPUTEKLMODESMESH  Eigenmodes of a Gaussian kernel on mesh centroids.
%
%   [V, Vx, Vy, lambda] = computeKLmodesMesh(mesh, L)
%   [V, Vx, Vy, lambda] = computeKLmodesMesh(mesh, L, m)
%
%   Builds the Gaussian kernel k(r) = exp(-(r/L)^2) between all cell
%   centroids, takes its leading eigenvectors (Karhunen-Loeve modes) by
%   SVD, and evaluates the x and z derivatives of each mode by
%   differentiating the kernel expansion
%       v_k(x) = (1/lambda_k) * sum_j k(x, x_j) v_k(x_j).
%   Modes are returned in order of decreasing lambda.
%
% INPUTS
%   mesh    object or struct with field
%           .xc      [M x 2] cell centroids [x, z] (m)
%   L       scalar kernel length scale (m)
%
% OPTIONS
%   m       number of modes to return; empty or omitted for all  [all]
%
% OUTPUTS
%   V       [M x m] modes, orthonormal columns
%   Vx      [M x m] derivative of each mode with respect to x (1/m)
%   Vy      [M x m] derivative of each mode with respect to the second
%           coordinate, z (1/m)
%   lambda  [m x 1] singular values of the kernel (the eigenvalues, since
%           the kernel is symmetric positive semi-definite)
%
% NOTES
%   - The full set uses svd(K, 'econ'); a truncated set uses svds. K is a
%     dense M x M matrix, so memory grows as M^2.
%   - Modes with lambda = 0 get zero derivatives.
%   - Not called by the driver scripts.
%
% SEE ALSO
%   svd, svds
%
% AUTHORS
%   Rishav Mallick, 2025

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
