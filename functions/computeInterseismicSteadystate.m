function [sol_f,sol_22,sol_23] = computeInterseismicSteadystate(rcv,shz,Vpl,lambda)
% Earthquake cycle calculation using semi-analytical solutions to elastic
% boundary value problem
% INPUTS:
% rcv,shz - geometry objects for fault & shear zone
% Vpl - relative plate convergence velocity
% lambda - regularization parameter (typically in the range 0.1 to 10,
%                                    default value if 10)
% OUTPUTS:
% sol_f - BEM solution for fault
% sol_22,sol_23 - BEM solution for two deviatoric components (viscous shearzones)
% 
% AUTHOR:
% Rishav Mallick, JPL 2024

rcv.Vpl = Vpl.*ones(rcv.N,1);

% Long-term strain rate calculation
[e22_dev, e23] = getStrainratesLongterm(shz,rcv.dip(1)*pi/180,[0,20e3],[-140e3,35e3]);
shz.e22pl = e22_dev.*Vpl;% 1/s
shz.e23pl = -e23.*Vpl;% 1/s
%% compute stress interaction kernels
% evl contains the following as N-d matrices
% KK - fault-fault interactions [rcv.N x rcv.N]
% KL - fault-shz interactions [shz.N x rcv.N x 2]
% LK - shz-fault interactions [rcv.N x shz.N x 2]
% LL - shz-shz interactions [shz.N x shz.N x 2 x 2]

% evl_orig = computeAllStressKernelsBem(rcv,shz,boundary,'kernelmodify',0);
% evl = computeAllStressKernelsBem(rcv,shz,boundary,'kernelmodify',1);
load('kernels/evl_orig.mat','evl_orig');

%% compute quasi-steady rates for 

% define locked zone on megathrust
locked = rcv.pinnedPosition;

fullstresskernel = [evl_orig.KK,         evl_orig.LK(:,:,1),         evl_orig.LK(:,:,2);...
                evl_orig.KL(:,:,1),      evl_orig.LL(:,:,1,1),       evl_orig.LL(:,:,1,2);...
                evl_orig.KL(:,:,2),      evl_orig.LL(:,:,2,1),       evl_orig.LL(:,:,2,2)];
% construct RHS (long-term loading rate due to locked section
bfit = fullstresskernel*[rcv.Vpl;shz.e22pl;shz.e23pl];
bfit(1:length(find(locked))) = [];% need to remove locked part of the domain


stresskernel = [evl_orig.KK(~locked,~locked),  evl_orig.LK(~locked,:,1),   evl_orig.LK(~locked,:,2);...
                evl_orig.KL(:,~locked,1),      evl_orig.LL(:,:,1,1),       evl_orig.LL(:,:,1,2);...
                evl_orig.KL(:,~locked,2),      evl_orig.LL(:,:,2,1),       evl_orig.LL(:,:,2,2)];


%% compute regularized BEM solution
% lambda = 10;

regularize_matrix = eye(shz.N*2+rcv.N-length(find(locked)));
regularize_matrix(1:length(find(~locked)),:) = [];

% solve BEM equations
sol_tot = ([stresskernel;lambda*regularize_matrix])\([bfit;zeros(shz.N*2,1)]);

% extract fault slip rate
sol_f = zeros(rcv.N,1);
sol_f(~locked) = sol_tot(1:length(find(~locked)));
% extract shear zone strain rates
sol_v = sol_tot(length(find(~locked))+1:end);
sol_22 = sol_v(1:end/2);
sol_23 = sol_v(end/2+1:end);

end

