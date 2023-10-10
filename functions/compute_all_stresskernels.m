function evl = compute_all_stresskernels(rcv,shz,boundary)
% Function that takes in a given shear zone data structure 'shz', and fault 'rcv', computes all relevant stress kernels
% 
% INPUTS
% rcv - object or data structure that contains fault geometry and Elastic parameters
% shz - object or data structure that contains shear zone geometry and Elastic parameters
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

mu = rcv.earthModel.G;
nu = rcv.earthModel.nu;

% create data structure to hold all kernels
evl = [];

%% compute fault-fault interaction kernels
disp('Computing fault - fault traction kernels')

% [K,~] = computetractionkernels(rcv,rcv);
[K,~] = computeFaultTractionKernelsBem(rcv,rcv,boundary);

% store shear traction kernel in 'evl'
evl.KK = K;

%% compute fault-shz deviatoric interaction kernel
% initialize stress kernels
LL1 = zeros(shz.N,rcv.N);
LL2 = zeros(shz.N,rcv.N);
LL3 = zeros(shz.N,rcv.N);

disp('Computing fault - shear zone traction kernels')
% xc = shz.xc;% here we don't need to flip the sign of z-coordinate
% for i = 1:rcv.N
%     m = [rcv.x(i,2) rcv.x(i,1) rcv.W(i) rcv.dip(i) 1];
%     [s22,s23,s33] = EdgeStress(m,xc(:,1),xc(:,2),nu,mu);
%     LL1(:,i) = s22;
%     LL2(:,i) = s23;
%     LL3(:,i) = s33;
% end

[LL1,LL3,LL2] = computeFaultTractionKernelsBem(rcv,shz,boundary);

% initialize deviatoric K-L kernel [Nshz x Nf x 2]
evl.KL = zeros(shz.N,rcv.N,2);
evl.KL(:,:,1) = (LL1-LL3)./2;
evl.KL(:,:,2) = LL2;

%% compute shz-shz interaction kernels
% initialize stress kernels
LL1 = zeros(shz.N,shz.N);
LL2 = zeros(shz.N,shz.N);
LL3 = zeros(shz.N,shz.N);

% this are 3x3 stress kernels
LL = zeros(shz.N,shz.N,3,3);

% deviatoric stress kernels [2x2]
evl.LL = zeros(shz.N,shz.N,2,2);

% source strain 100;010;001
I = eye(3);

tic
disp('Computing shear zone - shear zone stress kernels')
% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);
A = shz.A;A(:,2) = -A(:,2);
B = shz.B;B(:,2) = -B(:,2);
C = shz.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k = 1:shz.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);
        
        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        LL1(:,k) = s22(:);
        LL2(:,k) = s23(:);
        LL3(:,k) = s33(:);

        % LL(:,k,1,i)=s22(:);
        % LL(:,k,2,i)=s23(:);
        % LL(:,k,3,i)=s33(:);
    end
    LL(:,:,1,i) = LL1;
    LL(:,:,2,i) = LL2;
    LL(:,:,3,i) = LL3;
end

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
Evals_corrected(real(Evals_corrected)>0) = 0;

% reconstruct deviatoric stress kernel
L_deviatoric_corrected = real(Evectors*diag(Evals_corrected)/Evectors);

% store deviatoric kernels in 'evl'
evl.LL(:,:,1,1) = L_deviatoric_corrected(1:end/2,1:end/2);
evl.LL(:,:,1,2) = L_deviatoric_corrected(1:end/2,end/2+1:end);
evl.LL(:,:,2,1) = L_deviatoric_corrected(end/2+1:end,1:end/2);
evl.LL(:,:,2,2) = L_deviatoric_corrected(end/2+1:end,end/2+1:end);
toc

%% compute stress interactions from shear zones to fault (LK kernels)

% evaluate stress at fault center
xc = rcv.xc; 
xc(:,2) = -xc(:,2);

% only need 1 component of stress for faults
LL1 = zeros(rcv.N,shz.N);

% full stress kernels [N_fault x N_shz x 3]
LL = zeros(rcv.N,shz.N,3);

tic
disp('Computing shear zone - fault traction kernels')
for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k = 1:shz.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);
        % compute traction vector for fault plane orientation
        t=[s22.*rcv.nv(:,1)+s23.*rcv.nv(:,2), ...
            s23.*rcv.nv(:,1)+s33.*rcv.nv(:,2)];
        % rotate traction vector to fault-shear direction
        LL1(:,k) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);        

    end
    LL(:,:,i) = LL1;
end

% due to deviatoric state, e22 = -e33. Incorporating this into the kernels
% shape of kernel: [Nf x Nshz x 2]
evl.LK = zeros(rcv.N,shz.N,2);
evl.LK(:,:,1) = LL(:,:,1) - LL(:,:,3);
evl.LK(:,:,2) = LL(:,:,2);
toc

end


















