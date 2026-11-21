function evl = certifyShzKernel(evl, shz, varargin)
% CERTIFYSHZKERNEL  Make the mantle-mantle kernel provably free of growing modes.
%
%   evl = certifyShzKernel(evl, shz)
%   evl = certifyShzKernel(evl, shz, 'Tolerance', tol)
%
%   For a cell-averaged (Galerkin) kernel, the area-weighted deviatoric
%   mantle-mantle kernel S = diag([a; a]) * [LL11 LL12; LL21 LL22] is,
%   exactly, symmetric (Betti reciprocity) and negative semidefinite (the
%   elastic energy of any eigenstrain distribution is non-negative). With
%   quadrature the symmetry holds to ~1%. This function
%     1. symmetrises S (removes only quadrature error),
%     2. checks that every eigenvalue of the symmetric S is negative, and
%        raises an error otherwise,
%   and rebuilds LL = diag([a; a]) \ S. Nothing else is changed.
%
%   Guarantee: with the fault locked, every mantle mode of the certified
%   kernel decays, for any viscosity distribution and any power-law
%   exponent at any stress state (no growing and no neutral modes). The
%   linearised mantle system is L*D with D block-diagonal and positive
%   definite in each cell, and L*D is similar to X*S*X with X symmetric
%   positive definite, whose eigenvalues are real and < 0.
%
% INPUTS
%   evl     kernel struct from computeAllStressKernelsBem (uses .LL)
%   shz     geometry.shearZoneReceiver object (uses .area)
%
% OPTIONS
%   'Tolerance'   relative size of the largest eigenvalue, above which the
%                 kernel is called indefinite rather than having
%                 zero-energy modes (only the error message differs)  [1e-3]
%                 From the measured margin: zero-energy modes perturbed by
%                 3-point quadrature reach 1.4e-4 (criss-cross mesh), while
%                 the centroid example kernel has 4.3e-2.
%
% OUTPUTS
%   evl     same struct with .LL certified and .certificate:
%           .asymmetry       ||S - S'|| / ||S|| before symmetrisation
%           .maxEigenvalue   largest eigenvalue of the symmetric S, divided
%                            by its largest magnitude (negative)
%           .tolerance       tolerance used
%
% NOTES
%   - Centroid-collocated kernels ('cellaverage', 1) are not symmetric or
%     negative semidefinite and fail the check (the example mesh has 65
%     positive eigenvalues in the symmetric part); this is intended.
%   - Meshes that admit zero-energy modes (compatible piecewise-constant
%     strains that cause no stress, e.g. criss-cross meshes) fail with
%     'certifyShzKernel:zeroEnergyModes': quadrature moves those zero
%     eigenvalues slightly off zero, and such modes neither decay nor
%     grow, so the periodic cycle is not determined. Use another mesh.
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
    error('certifyShzKernel:notNegative', ...
        ['The area-weighted mantle-mantle kernel is not negative semidefinite: ' ...
         '%d eigenvalues of its symmetric part are positive, the largest %.2e of the ' ...
         'spectrum (tolerance %.0e). Use cell-averaged kernels (''cellaverage'', 3).'], ...
        sum(lambda > tol*scale), maxEig, tol);
end

if maxEig >= 0
    error('certifyShzKernel:zeroEnergyModes', ...
        ['The area-weighted mantle-mantle kernel has %d eigenvalues of its symmetric part at or ' ...
         'just above zero (largest %.2e of the spectrum, below the indefinite threshold %.0e): ' ...
         'zero-energy modes of the mesh (e.g. a criss-cross mesh), moved off zero by quadrature. ' ...
         'They neither decay nor grow, so the cycle is not determined. Use a mesh without them.'], ...
        sum(lambda >= 0), maxEig, tol);
end

L = S ./ w;
evl.LL(:,:,1,1) = L(1:M, 1:M);
evl.LL(:,:,1,2) = L(1:M, M+1:end);
evl.LL(:,:,2,1) = L(M+1:end, 1:M);
evl.LL(:,:,2,2) = L(M+1:end, M+1:end);
evl.certificate = struct('asymmetry', asymmetry, 'maxEigenvalue', maxEig, ...
    'tolerance', tol);
end
