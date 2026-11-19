function evl = computeAllStressKernelsBem(rcv,shz,boundary,varargin)
% COMPUTEALLSTRESSKERNELSBEM  Assemble all stress interaction kernels.
%
%   evl = computeAllStressKernelsBem(rcv, shz, boundary)
%   evl = computeAllStressKernelsBem(rcv, shz, boundary, 'kernelmodify', 0)
%
%   Computes the four blocks of the elastic interaction operator between
%   fault patches (rcv) and mantle shear-zone cells (shz), corrected for
%   the model boundary by the boundary element method: fault slip to fault
%   shear traction (KK), fault slip to mantle deviatoric stress (KL),
%   mantle deviatoric strain to fault shear traction (LK) and mantle
%   deviatoric strain to mantle deviatoric stress (LL). Mantle strain is
%   deviatoric (e33 = -e22), so each cell has two strain components
%   (e22, e23) and two stress components, sigma22' = (sxx - szz)/2 and
%   sigma23 = sxz.
%
% INPUTS
%   rcv       geometry.receiver object: megathrust patches (N patches)
%   shz       geometry.shearZoneReceiver object: mantle cells (M cells)
%   boundary  geometry.receiver object: model boundary patches
%             .Vpl  0: zero displacement; 1: zero traction (.Vx, .Vz
%             are ignored: the kernels solve the homogeneous problem)
%
% OPTIONS
%   kernelmodify  1: set the eigenvalues of the deviatoric LL kernel with
%                 positive real part to 0; 0: keep LL as computed  [1]
%
% OUTPUTS
%   evl       struct with fields
%             .KK  [N x N] shear traction on rcv per unit slip on rcv
%                  (MPa/m)
%             .KL  [M x N x 2] deviatoric stress in shz per unit slip on
%                  rcv; (:,:,1) = (sxx - szz)/2, (:,:,2) = sxz (MPa/m)
%             .LK  [N x M x 2] shear traction on rcv per unit deviatoric
%                  strain in shz; (:,:,1) from e22 (with e33 = -e22),
%                  (:,:,2) from e23 (MPa)
%             .LL  [M x M x 2 x 2] deviatoric stress in shz per unit
%                  deviatoric strain in shz; 3rd index = receiver stress
%                  component (sigma22', sigma23), 4th index = source strain
%                  component (e22, e23) (MPa)
%
% NOTES
%   - Kernel modification (default): the 2M x 2M matrix
%     [LL11 LL12; LL21 LL22] is eigendecomposed, every eigenvalue with
%     positive real part is set to 0, and the matrix is rebuilt (real part
%     kept). Positive eigenvalues come from the discretisation and would
%     make the viscoelastic system grow instead of relax. Call with
%     'kernelmodify', 0 for the unmodified kernels (evl_orig).
%   - KK, KL and LK are never modified.
%   - A comment in the code says the index construction of the deviatoric
%     LL kernel still needs to be checked.
%   - Units assume the shear modulus of the earth model is in MPa.
%   - Prints progress and timings to the command window.
%
% SEE ALSO
%   computeFaultTractionKernelsBem, computeShzStressKernelsBem,
%   computeAllDisplacementKernelsBem, assembleStressKernel,
%   negateGrowingModes
%
% AUTHORS
%   Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

% Check the number of optional arguments
numOptionalArgs = nargin - 3; % Subtract the three required arguments

% default is to modify kernels
modifykernel = 1;
% Parse name-value pairs
for i = 1:2:numOptionalArgs
    argname = varargin{i};
    argvalue = varargin{i + 1};

    % Check the name and assign the value accordingly
    switch lower(argname)
        case 'kernelmodify'
            modifykernel = argvalue;        
        otherwise
            error(['Unknown parameter: ', argname]);
    end
end

% create data structure to hold all kernels
evl = [];

%% compute fault-fault interaction kernels
disp('Computing fault - fault traction kernels')

[K,~] = computeFaultTractionKernelsBem(rcv,rcv,boundary);

% store shear traction kernel in 'evl'
evl.KK = K;

%% compute fault-shz deviatoric interaction kernel

disp('Computing fault - shear zone traction kernels')
% note the ordering of kernels here is [sxx,szz,sxz]
[LL1,LL3,LL2] = computeFaultTractionKernelsBem(rcv,shz,boundary);

% initialize deviatoric K-L kernel [Nshz x Nf x 2]
evl.KL = zeros(shz.N,rcv.N,2);
evl.KL(:,:,1) = (LL1-LL3)./2;
evl.KL(:,:,2) = LL2;

%% compute shz-shz interaction kernels

disp('Computing shear zone - shear zone stress kernels')
tic
LL = computeShzStressKernelsBem(shz,shz,boundary);

% deviatoric stress kernels [2x2]
evl.LL = zeros(shz.N,shz.N,2,2);

% construct deviatoric stress kernel and remove positive eigen values
% here the first 2 indices are for eigen strain source
%           last 2 indices are for the stress component
% evl.LL2333 corresponds to e23 source resulting in a change in s33
% in the large 3x3 matrix this would be in position [3-row,2-column]
% L2222 = (evl.LL2222 - evl.LL3322 - evl.LL2233 + evl.LL3333)./2;
% L2322 = (evl.LL2322 - evl.LL2333)./2;
% L2223 = (evl.LL2223 - evl.LL3323);
% L2323 = evl.LL2323;

% need to check/verify this construction (the indices always confuse me)
L11 = (LL(:,:,1,1) - LL(:,:,1,3) - LL(:,:,3,1) + LL(:,:,3,3))./2;
L12 = (LL(:,:,1,2) - LL(:,:,3,2))./2;
L21 = (LL(:,:,2,1) - LL(:,:,2,3));
L22 = LL(:,:,2,2);

L_deviatoric = [L11 L12;...
        L21 L22];
[Evectors,Evals] = eig(L_deviatoric);
% remove eigen values that cause instabilities 
Evals_corrected = diag(Evals);

if modifykernel == 1
    Evals_corrected(real(Evals_corrected)>0) = 0;
end

% reconstruct deviatoric stress kernel
L_deviatoric_corrected = real(Evectors*diag(Evals_corrected)/Evectors);

% store deviatoric kernels in 'evl'
evl.LL(:,:,1,1) = L_deviatoric_corrected(1:end/2,1:end/2);
evl.LL(:,:,1,2) = L_deviatoric_corrected(1:end/2,end/2+1:end);
evl.LL(:,:,2,1) = L_deviatoric_corrected(end/2+1:end,1:end/2);
evl.LL(:,:,2,2) = L_deviatoric_corrected(end/2+1:end,end/2+1:end);
toc

%% compute stress interactions from shear zones to fault (LK kernels)

disp('Computing shear zone - fault traction kernels')
tic
LL = computeShzStressKernelsBem(shz,rcv,boundary);

% due to deviatoric state, e22 = -e33. Incorporating this into the kernels
% shape of kernel: [Nf x Nshz x 2]
evl.LK = zeros(rcv.N,shz.N,2);
% evl.LK(:,:,1) = LL(:,:,1) - LL(:,:,3);
% evl.LK(:,:,2) = LL(:,:,2);

evl.LK(:,:,1) = LL(:,:,1,1) - LL(:,:,1,3);
evl.LK(:,:,2) = LL(:,:,1,2);
toc

end


















