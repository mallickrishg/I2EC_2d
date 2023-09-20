function devl = compute_all_dispkernels(obs,rcv,shz)
% Function that takes in a given shear zone data structure 'shz', and fault 'rcv', computes all relevant displacement kernels
% 
% INPUTS
% obs - array containing the x2,x3 coordinates of the observation points 
% rcv - object or data structure that contains fault geometry and Elastic parameters
% shz - object or data structure that contains shear zone geometry and Elastic parameters
% 
% OUTPUTS
% devl - data structure with the following kernels
% KO  - displacement kernels due to fault slip 
% KO(:,:,1) - horizontal component of displacement 
% KO(:,:,2) - vertical component of displacement
% LO  - displacement kernels due to strain in shear zones 
% LO(:,:,1,1) and LO(:,:,1,2) - horizontal component of displacement 
% LO(:,:,2,1) and LO(:,:,2,2) - vertical  component of displacement 
%
% Authors:
% Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

mu = rcv.earthModel.G;
nu = rcv.earthModel.nu;

% create data structure to hold disp kernels
devl = [];

%% compute fault-fault interaction kernels
disp('Computing displacement kernels')

Nobs=length(obs(:,1));
%% compute fault displacement kernels

Gf_d = computedisplacementGFs(rcv,obs);

% Nobs x shz.N x 2 corresponding to the horizontal and 
KO = zeros(Nobs,shz.N,2);
devl.KO(:,:,1) = Gf_d(1:2:end,:);
devl.KO(:,:,2) = Gf_d(2:2:end,:);

%% Displacement kernels for shear zones 

% initialize displacement kernels
LO1 = zeros(Nobs,shz.N);
LO2 = zeros(Nobs,shz.N);

% 2x3 displacemet kernels: 3 source strains and 2 components of
% displacement 
LO = zeros(Nobs,shz.N,2,3);

% displacement kernels after incorporating deviatoric stress constraints [2x2]
devl.LO = zeros(Nobs,shz.N,2,2);

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
        LO2(:,k)= u3(:);
      
    end
    LO(:,:,1,i) = LO1; % horizontal component of displacement 
    LO(:,:,2,i) = LO2; % vertical component of displacement
end

% due to deviatoric state, e22 = -e33. Incorporating this into the kernels
% shape of kernel: [Nobs x Nshz x 2 x 2]
devl.LO(:,:,1,1) = LO(:,:,1,1) - LO(:,:,1,3); 
devl.LO(:,:,2,1) = LO(:,:,2,1) - LO(:,:,2,3);
devl.LO(:,:,1,2)=  LO(:,:,1,2);
devl.LO(:,:,2,2)=  LO(:,:,2,2);

%%



