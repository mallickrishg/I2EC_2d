function evl = compute_all_stresskernels(rcv,shz)
% for a given shear zone data structure 'shz', and fault 'rcv', compute all
% relevant stress kernels
% Rishav Mallick & Sharadha Sathiakumar, 2023

mu = rcv.earthModel.G;
nu = rcv.earthModel.nu;

% create data structure to hold all kernels
evl = [];

%% compute fault-fault interaction kernels
[K,~] = computetractionkernels(rcv,rcv);
evl.KK = K;

%% compute shz-shz interaction kernels
% initialize stress kernels
LL1 = zeros(shz.N,shz.N);
LL2 = zeros(shz.N,shz.N);
LL3 = zeros(shz.N,shz.N);
LL = zeros(shz.N,shz.N,3,3);

% deviatoric stress kernels
evl.LL = zeros(shz.N,shz.N,2,2);

% source strain 100;010;001
I=eye(3);
tic

% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);
A = shz.A;A(:,2) = -A(:,2);
B = shz.B;B(:,2) = -B(:,2);
C = shz.C;C(:,2) = -C(:,2);

for i=1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k=1:shz.N
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

toc

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

evl.LL(:,:,1,1) = L_deviatoric_corrected(1:end/2,1:end/2);
evl.LL(:,:,1,2) = L_deviatoric_corrected(1:end/2,end/2+1:end);
evl.LL(:,:,2,1) = L_deviatoric_corrected(end/2+1:end,1:end/2);
evl.LL(:,:,2,2) = L_deviatoric_corrected(end/2+1:end,end/2+1:end);

end


















