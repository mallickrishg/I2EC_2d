function stresskernel = assembleStressKernel(evl,locked)
% ASSEMBLESTRESSKERNEL  Assemble the stress kernel of the linear ODE system.
%
%   stresskernel = assembleStressKernel(evl, locked)
%
%   Stacks the fault and shear-zone interaction kernels into one square
%   matrix that maps [slip rate on creeping patches; e22 rate; e23 rate]
%   to [shear stressing rate on creeping patches; s22 rate; s23 rate].
%   Locked patches are removed from rows and columns. Dividing the rows by
%   the viscosities gives the operator A of dy/dt = A*y + c used by
%   solveAnalyticIvp:
%       A = diag(1./viscosityvector)*stresskernel
%
% INPUTS
%   evl     struct of stress kernels from computeAllStressKernelsBem
%           .KK      [N x N] fault-to-fault shear traction (MPa/m)
%           .KL      [M x N x 2] fault-to-mantle deviatoric stress (MPa/m)
%           .LK      [N x M x 2] mantle-to-fault shear traction (MPa per
%                    unit strain)
%           .LL      [M x M x 2 x 2] mantle-to-mantle deviatoric stress
%                    (MPa per unit strain)
%   locked  [N x 1] logical, true for locked fault patches
%
% OUTPUTS
%   stresskernel  [Nc+2M x Nc+2M] assembled kernel, where Nc is the number
%                 of creeping (unlocked) patches
%
% NOTES
%   - The drivers pass the default kernels of computeAllStressKernelsBem
%     (cell-averaged and certified, saved as examples/kernels/evl_orig.h5), which
%     have no growing modes; negateGrowingModes in the solvers then changes
%     nothing. With centroid kernels it negates their growing modes.
%
% SEE ALSO
%   computeAllStressKernelsBem, negateGrowingModes, solveAnalyticIvp
%
% AUTHORS
%   Rishav Mallick, 2025

stresskernel = [evl.KK(~locked,~locked),  evl.LK(~locked,:,1),   evl.LK(~locked,:,2);...
                evl.KL(:,~locked,1),      evl.LL(:,:,1,1),       evl.LL(:,:,1,2);...
                evl.KL(:,~locked,2),      evl.LL(:,:,2,1),       evl.LL(:,:,2,2)];

end