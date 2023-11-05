function LL = computeShzStressKernels(src,shz)
% half-space stress kernel computation in 2-d plane strain returns stress/traction kernels: 
% each kernel is a component of the full 2-d stress tensor [sxx,sxz,szz] 
%                                 or fault-shear tractions [tau_s]
% in response to each non-elastic strain component [exx,exz,ezz] 
% (note: this ordering of components is different from  computeFaultTractionKernels)
% 
% OUTPUTS:
% LL is a [Nshz x Nsrc x 3 x 3] matrix (when shz is shearZoneReceiver)
% LL is a [Nshz x Nsrc x 3] matrix (when shz is receiver)
% for shearZoneReceiver: 
%          4th index is for eigen strain source
%          3rd index is for the stress component
% if we want a kernel for - 
% exz source resulting in szz
% in the N x M X 3 x 3 matrix this would be
% position [N x M x 3-row,2-column]
% for receiver:
%          4th index is for eigen strain source
%          3rd index is for the traction component
% if we want a kernel for - 
% ezz source resulting in tau_n
% in the N x M X 2 x 3 matrix this would be 
% position [N x M x 2-row,3-column]
% 
% Author:
% Rishav Mallick, JPL, 2023

mu = src.earthModel.G;
nu = src.earthModel.nu;

% initialize stress kernels
LL1 = zeros(shz.N,src.N);
LL2 = zeros(shz.N,src.N);
LL3 = zeros(shz.N,src.N);

% this are 3x3 stress kernels
LL = zeros(shz.N,src.N,3,3);

% source strain 100;010;001
I = eye(3);

% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);

A = src.A;A(:,2) = -A(:,2);
B = src.B;B(:,2) = -B(:,2);
C = src.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to exx source
    % i = 2 corresponds to exz source
    % i = 3 corresponds to ezz source
    parfor k = 1:src.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);
        
        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        if isa(shz,'geometry.shearZoneReceiver')
            LL1(:,k) = s22(:);
            LL2(:,k) = s23(:);
            LL3(:,k) = s33(:);
        elseif isa(shz,'geometry.receiver')
            % compute traction vector for fault plane orientation
            t=[s22.*shz.nv(:,1) - s23.*shz.nv(:,2), ...
                -s23.*shz.nv(:,1) + s33.*shz.nv(:,2)];
            % rotate traction vector to fault-shear & fault-normal direction
            LL1(:,k) = shz.dv(:,1).*t(:,1) + shz.dv(:,2).*t(:,2);
            LL2(:,k) = shz.nv(:,1).*t(:,1) + shz.nv(:,2).*t(:,2);
        else
            LL1(:,k) = s22(:);
            LL2(:,k) = s23(:);
            LL3(:,k) = s33(:);
            % error('not a recognized geometry. provide either geometry.shearZoneReceiver or geometry.receiver')
        end

    end
    if isa(shz,'geometry.shearZoneReceiver')
        LL(:,:,1,i) = LL1; % sxx
        LL(:,:,2,i) = LL2; % sxz
        LL(:,:,3,i) = LL3; % szz
    elseif isa(shz,'geometry.receiver')
        LL(:,:,1,i) = LL1; % tau_shear
        LL(:,:,2,i) = LL2; % tau_normal
    else
        LL(:,:,1,i) = LL1; % sxx
        LL(:,:,2,i) = LL2; % sxz
        LL(:,:,3,i) = LL3; % szz
        % error('not a recognized geometry. provide either geometry.shearZoneReceiver or geometry.receiver')
    end
end

end