function LL = computeShzStressKernels(src,shz)
% half-space stress kernel computation in 2-d plane strain returns all 9 stress kernels: 
% each kernel is a component of the full 2-d stress tensor [sxx,sxz,szz]
% in response to each strain component [exx,exz,ezz] 
% (note: this ordering is different from the fault kernels)
% 
% OUTPUTS:
% LL is a [Nshz x Nsrc x 3 x 3] matrix
% here the 4th index is for eigen strain source
%          3rd index is for the stress component
% if we want a kernel for - 
% exz source resulting in szz
% in the N x M X 3 x 3 matrix this would be
% position [N x M x 3-row,2-column]
% 
% Author:
% Rishav Mallick, JPL, 2023

% initialize stress kernels
LL1 = zeros(shz.N,src.N);
LL2 = zeros(shz.N,src.N);
LL3 = zeros(shz.N,src.N);

% this are 3x3 stress kernels
LL = zeros(shz.N,shz.N,3,3);

% source strain 100;010;001
I = eye(3);


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
    parfor k = 1:src.N
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

    end
    LL(:,:,1,i) = LL1;
    LL(:,:,2,i) = LL2;
    LL(:,:,3,i) = LL3;
end

end