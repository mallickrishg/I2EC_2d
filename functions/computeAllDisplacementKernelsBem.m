function devl = computeAllDisplacementKernelsBem(obs,rcv,shz,boundary,scalar_value)
% COMPUTEALLDISPLACEMENTKERNELSBEM  Assemble displacement kernels at observation points.
%
%   devl = computeAllDisplacementKernelsBem(obs, rcv, shz, boundary, scalar_value)
%
%   Computes the displacement at observation points due to slip on the
%   fault patches (KO) and deviatoric strain in the mantle cells (LO),
%   corrected for the model boundary by the boundary element method.
%
% INPUTS
%   obs           [Nobs x 2] observation points [x, z] (m)
%   rcv           geometry.receiver object: megathrust patches (N patches)
%   shz           geometry.shearZoneReceiver object: mantle cells (M cells)
%   boundary      geometry.receiver object: model boundary patches
%                 .Vpl  0: zero displacement; 1: zero traction
%   scalar_value  source amplitude of each kernel column (slip in m or
%                 strain); 1 gives kernels per unit source
%
% OUTPUTS
%   devl          struct with fields
%                 .KO  [Nobs x N x 2] displacement due to slip on each
%                      patch; (:,:,1) = ux, (:,:,2) = uz
%                 .LO  [Nobs x M x 2 x 2] displacement due to deviatoric
%                      strain in each cell; 3rd index = displacement
%                      component (ux, uz), 4th index = source strain
%                      component (e22 with e33 = -e22, e23)
%
% NOTES
%   - The driver scripts call it with scalar_value = 1.
%   - Prints progress to the command window.
%
% SEE ALSO
%   computeFaultDisplacementKernelsBem, computeShzDisplacementKernelsBem,
%   computeAllStressKernelsBem
%
% AUTHORS
%   Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

% create data structure to hold disp kernels
devl = [];

%% compute fault-fault interaction kernels
disp('Computing displacement kernels')

Nobs=length(obs(:,1));
%% compute fault displacement kernels

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



