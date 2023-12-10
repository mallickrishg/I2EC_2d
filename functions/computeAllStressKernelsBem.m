function evl = computeAllStressKernelsBem(rcv,shz,boundary,varargin)
% Function that takes in a given shear zone data structure 'shz', and fault 'rcv', computes all relevant stress kernels
% 
% INPUTS
% rcv - object or data structure that contains fault geometry and Elastic parameters
% shz - object or data structure that contains shear zone geometry and Elastic parameters
% boundary - object or data structure that contains the boundary mesh and Elastic parameters
% 
% OUTPUTS
% evl - data structure with the following kernels
% KK  - rcv(source)-shz(receiver) traction kernel (only in fault-shear direction) [Nf x Nf]
% KL  - rcv(source)-shz(receiver) deviatoric stress kernel [Nshz x Nf x 2]
% LK  - shz(source)-rcv(receiver) traction kernel (sources are deviatoric
%       eigen strains, receiver traction is only in shear component) [Nf x Nshz x 2]
% LL  - shz(source)-shz(receiver) deviatoric strain source and 
%       deviatoric receiver stress [Nshz x Nshz x 2 x 2]
% 
% Authors:
% Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

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


















