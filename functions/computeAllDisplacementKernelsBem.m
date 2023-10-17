function devl = computeAllDisplacementKernelsBem(obs,rcv,shz,boundary,scalar_value)
% Function that takes in a given shear zone data structure 'shz', and fault 'rcv', computes all relevant displacement kernels
% 
% INPUTS
% obs - array containing the x2,x3 coordinates of the observation points 
% rcv - object or data structure that contains fault geometry and Elastic parameters
% shz - object or data structure that contains shear zone geometry and Elastic parameters
% scalar_value - rescale BEM problem (typically the plate velocity)
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

% create data structure to hold disp kernels
devl = [];

%% compute fault-fault interaction kernels
disp('Computing displacement kernels')

Nobs=length(obs(:,1));
%% compute fault displacement kernels

% Gf_d = computedisplacementGFs(rcv,obs);
[Gx_d,Gz_d] = computeFaultDisplacementKernelsBem(rcv,obs,boundary,scalar_value);

% Nobs x shz.N x 2 corresponding to the horizontal and vertical kernels
devl.KO = zeros(Nobs,rcv.N,2);
devl.KO(:,:,1) = Gx_d;
devl.KO(:,:,2) = Gz_d;

%% Displacement kernels for shear zones 

% 2x3 displacemet kernels: 
% 3 source strains and 2 components of displacement 

% displacement kernels after incorporating deviatoric stress constraints [2x2]
devl.LO = zeros(Nobs,shz.N,2,2);

LO = computeShzDisplacementKernelsBem(shz,obs,boundary,scalar_value);

% due to deviatoric state, e22 = -e33. Incorporating this into the kernels
% shape of kernel: [Nobs x Nshz x 2 x 2]
devl.LO(:,:,1,1) = LO(:,:,1,1) - LO(:,:,1,3); % ux
devl.LO(:,:,2,1) = LO(:,:,2,1) - LO(:,:,2,3); % uz
devl.LO(:,:,1,2)=  LO(:,:,1,2); % ux
devl.LO(:,:,2,2)=  LO(:,:,2,2); % uz

%%



