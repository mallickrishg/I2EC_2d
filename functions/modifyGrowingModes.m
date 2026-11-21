function evl = modifyGrowingModes(evl, mode)
% MODIFYGROWINGMODES  Remove the growing modes of the mantle-mantle kernel.
%
%   evl = modifyGrowingModes(evl, mode)
%
%   Eigendecomposes the 2M x 2M deviatoric mantle kernel
%   [LL11 LL12; LL21 LL22], changes every eigenvalue with positive real
%   part, and rebuilds LL (real part kept). The collocated kernel has a few
%   such eigenvalues on every mesh tried; they would make the viscous
%   relaxation grow instead of decay.
%
%   The earlier, empirical fix for centroid kernels ('cellaverage', 1), for
%   numerical (ODE) solves only. The default kernels of
%   computeAllStressKernelsBem are cell-averaged and certified and have no
%   growing modes, so they need no modification.
%
% INPUTS
%   evl     kernel struct from computeAllStressKernelsBem; only .LL is used
%   mode    0  keep LL as it is
%           1  'zero': set the growing eigenvalues to 0. The modes become
%              neutral: nothing restores them, so they can accumulate over
%              many earthquake cycles
%           2  'damp': set them to -max|lambda| (largest eigenvalue
%              magnitude of the kernel). The modes decay at the fastest
%              rate the kernel has
%
% OUTPUTS
%   evl     same struct with .LL modified; KK, KL and LK are unchanged
%
% NOTES
%   - The number of growing modes of the linearised mantle system does
%     not depend on the rheology: for a reciprocal kernel and any power-law
%     exponent and stress state it equals the number of positive
%     eigenvalues of LL. Modifying the kernel therefore works for n = 1
%     and n = 3 alike.
%   - Neither change is small at the cell scale; results for single cells
%     or sharp cell-to-cell contrasts should not be interpreted.
%
% EXAMPLE
%   evl = computeAllStressKernelsBem(rcv, shz, boundary, ...
%       'cellaverage', 1, 'kernelcertify', false);   % centroid kernel
%   evl = modifyGrowingModes(evl, 2);              % damped, for ode45
%
% SEE ALSO
%   computeAllStressKernelsBem, negateGrowingModes, odeViscoelastic

if ~(isscalar(mode) && ismember(mode, [0 1 2]))
    error('modifyGrowingModes:mode', ...
        'mode must be 0 (keep), 1 (zero) or 2 (damp), not %s.', mat2str(mode));
end
if mode == 0
    return
end

L = [evl.LL(:,:,1,1) evl.LL(:,:,1,2);
     evl.LL(:,:,2,1) evl.LL(:,:,2,2)];
[V, E] = eig(L);
lambda = diag(E);
growing = real(lambda) > 0;

if mode == 1
    lambda(growing) = 0;
else
    lambda(growing) = -max(abs(lambda));
end

L = real(V * diag(lambda) / V);

M = size(evl.LL, 1);
evl.LL(:,:,1,1) = L(1:M, 1:M);
evl.LL(:,:,1,2) = L(1:M, M+1:end);
evl.LL(:,:,2,1) = L(M+1:end, 1:M);
evl.LL(:,:,2,2) = L(M+1:end, M+1:end);
end
