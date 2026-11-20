function evl = certifyMantleKernel(evl, shz, varargin)
% CERTIFYMANTLEKERNEL  Make the mantle-mantle kernel provably free of growing modes.
%
%   evl = certifyMantleKernel(evl, shz)
%   evl = certifyMantleKernel(evl, shz, 'Tolerance', tol)
%
%   For a cell-averaged (Galerkin) kernel, the area-weighted deviatoric
%   mantle-mantle kernel S = diag([a; a]) * [LL11 LL12; LL21 LL22] is,
%   exactly, symmetric (Betti reciprocity) and negative semidefinite (the
%   elastic energy of any eigenstrain distribution is non-negative). With
%   quadrature the symmetry holds to ~1%. This function
%     1. symmetrises S (removes only quadrature error),
%     2. checks that the largest eigenvalue of the symmetric S is <= 0
%        within 'Tolerance' (default 1e-3) times its largest magnitude,
%        and raises an error otherwise,
%     3. sets eigenvalues that are positive only within the tolerance to
%        exactly 0 (zero-energy modes, e.g. compatible strain fields on
%        criss-cross meshes, perturbed by quadrature),
%   and rebuilds LL = diag([a; a]) \ S.
%
%   Guarantee: the certified kernel gives no growing modes for any
%   viscosity distribution and any power-law exponent at any stress
%   state. The linearised mantle system is L*D with D block-diagonal and
%   positive definite in each cell, and L*D is similar to X*S*X with X
%   symmetric positive definite, whose eigenvalues are real and <= 0.
%
% INPUTS
%   evl     kernel struct from computeAllStressKernelsBem (uses .LL)
%   shz     geometry.shearZoneReceiver object (uses .area)
%
% OPTIONS
%   'Tolerance'   relative tolerance on the largest eigenvalue  [1e-3]
%                 Chosen from the margin between quadrature noise and a
%                 truly indefinite kernel: zero-energy modes perturbed by
%                 3-point quadrature reach 1.4e-4 (criss-cross mesh), while
%                 the centroid example kernel has 4.3e-2. Clipping changes
%                 the spectrum by at most this fraction, below the
%                 quadrature error itself (antisymmetric part ~1e-2).
%
% OUTPUTS
%   evl     same struct with .LL certified and .certificate:
%           .asymmetry       ||S - S'|| / ||S|| before symmetrisation
%           .maxEigenvalue   largest eigenvalue of the symmetric S, divided
%                            by its largest magnitude, before clipping
%           .clipped         number of eigenvalues set from (0, tol] to 0
%           .tolerance       tolerance used
%
% NOTES
%   - Centroid-collocated kernels ('cellaverage', 1) are not symmetric or
%     negative semidefinite and fail the check (the example mesh has 65
%     positive eigenvalues in the symmetric part); this is intended.
%   - Only LL is changed. The fault blocks (KK, KL, LK) are not certified.
%
% SEE ALSO
%   computeAllStressKernelsBem

p = inputParser;
p.addParameter('Tolerance', 1e-3, @(t) isscalar(t) && t >= 0);
p.parse(varargin{:});
tol = p.Results.Tolerance;

M = size(evl.LL, 1);
w = [shz.area(:); shz.area(:)];
L = [evl.LL(:,:,1,1) evl.LL(:,:,1,2);
     evl.LL(:,:,2,1) evl.LL(:,:,2,2)];
S = w .* L;
asymmetry = norm(S - S', 'fro') / norm(S, 'fro');
S = (S + S') / 2;

lambda = eig(S);
scale = max(abs(lambda));
maxEig = max(lambda) / scale;
if maxEig > tol
    error('certifyMantleKernel:notNegative', ...
        ['The area-weighted mantle-mantle kernel is not negative semidefinite: ' ...
         '%d eigenvalues of its symmetric part are positive, the largest %.2e of the ' ...
         'spectrum (tolerance %.0e). Use cell-averaged kernels (''cellaverage'', 3).'], ...
        sum(lambda > tol*scale), maxEig, tol);
end

clipped = nnz(lambda > 0);
if clipped > 0
    [Q, E] = eig(S);
    e = diag(E);
    e(e > 0) = 0;
    S = Q * diag(e) * Q';
    S = (S + S') / 2;
end

L = S ./ w;
evl.LL(:,:,1,1) = L(1:M, 1:M);
evl.LL(:,:,1,2) = L(1:M, M+1:end);
evl.LL(:,:,2,1) = L(M+1:end, 1:M);
evl.LL(:,:,2,2) = L(M+1:end, M+1:end);
evl.certificate = struct('asymmetry', asymmetry, 'maxEigenvalue', maxEig, ...
    'clipped', clipped, 'tolerance', tol);
end
