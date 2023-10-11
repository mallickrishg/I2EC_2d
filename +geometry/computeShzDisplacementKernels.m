function LO = computeShzDisplacementKernels(shz,obs)
% Compute displacement kernels for shear zones 
% each kernel is a 2-component displacement field [ux,uz]
% in response to applied strain components [exx,exz,ezz]
% 
% OUTPUTS:
% Kdxx,Kdzz,Kdxz - for dip-slip source component
% Knxx,Knzz,Knxz - for tensile-slip source component
% 
% Author:
% Rishav Mallick, JPL, 2023
Nobs = length(obs(:,1));

% initialize displacement kernels
LO1 = zeros(Nobs,shz.N);
LO2 = zeros(Nobs,shz.N);

% 2 x 3 displacemet kernels: 3 source strains and 2 components of
% displacement 
LO = zeros(Nobs,shz.N,2,3);

% source strain 100;010;001
I = eye(3);

tic
disp('Computing disp kernels for shear zones')
% convert all depths to positive numbers
A = shz.A;A(:,2) = -A(:,2);
B = shz.B;B(:,2) = -B(:,2);
C = shz.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k = 1:shz.N
        [u2,u3]=computeDisplacementPlaneStrainTriangleShearZone( ...
            obs(:,1),obs(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            nu);

        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        LO1(:,k) = u2(:);
        LO2(:,k) = u3(:);
      
    end
    LO(:,:,1,i) = LO1; % horizontal component of displacement 
    LO(:,:,2,i) = LO2; % vertical component of displacement
end

end