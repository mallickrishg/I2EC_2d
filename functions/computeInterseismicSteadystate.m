function [sol_f,sol_22,sol_23] = computeInterseismicSteadystate(rcv,shz,evl_orig,Vpl,lambda)
% COMPUTEINTERSEISMICSTEADYSTATE  Steady interseismic rates with locking.
%
%   [sol_f, sol_22, sol_23] = computeInterseismicSteadystate(rcv, shz, ...
%                                                   evl_orig, Vpl, lambda)
%
%   Solves for the fully relaxed (t -> infinity) interseismic state of the
%   linear system when part of the megathrust is locked: the slip rate on
%   the creeping patches and the mantle strain rates for which the
%   stressing rate vanishes,
%       K(x_i; x_i)*[v; e22; e23] = K(x_i; x)*[Vpl; e22pl; e23pl],
%   where x_i are the unlocked degrees of freedom. K(x_i; x_i) is
%   singular, so the system is solved in the least-squares sense with a
%   zeroth-order Tikhonov penalty lambda*||[e22; e23]|| on the mantle
%   strain rates.
%
% INPUTS
%   rcv       geometry.receiver object: megathrust patches
%             .N               number of patches
%             .dip             [N x 1] dip (degrees); dip(1) is used
%             .pinnedPosition  [N x 1] logical, true for locked patches
%   shz       geometry.shearZoneReceiver object: mantle cells
%             .N               number of cells
%   evl_orig  struct of stress kernels (.KK, .KL, .LK, .LL) from
%             computeAllStressKernelsBem (default: certified, unmodified)
%   Vpl       scalar plate rate (m/s)
%   lambda    scalar Tikhonov weight, typically 0.1 to 10
%
% OUTPUTS
%   sol_f     [N x 1] fault slip rate (m/s); zero on locked patches
%   sol_22    [M x 1] deviatoric strain rate e22 (1/s)
%   sol_23    [M x 1] shear strain rate e23 (1/s)
%
% NOTES
%   - Stability: the regularisation removes the null space of K(x_i; x_i)
%     by pulling the mantle strain rates towards zero (not towards their
%     long-term values); fault slip rates are not penalised.
%   - The long-term strain rates come from getStrainratesLongterm with
%     fixed corner-flow offsets [0, 20e3] and [-140e3, 35e3] (m) and the
%     dip of the first patch.
%   - Side effects: sets rcv.Vpl, shz.e22pl and shz.e23pl on the (handle)
%     objects passed in.
%   - The locked rows are removed as the first numel(find(locked)) rows,
%     so the locked patches must be the first patches of rcv.
%   - Known issue: getStrainratesLongterm is called with four arguments,
%     but it needs a fifth (the power-law exponent n), so this function
%     errors as written.
%
% SEE ALSO
%   getStrainratesLongterm, computeAllStressKernelsBem, assembleStressKernel
%
% AUTHORS
%   Rishav Mallick (JPL), 2024

rcv.Vpl = Vpl.*ones(rcv.N,1);

% Long-term strain rate calculation
[e22_dev, e23] = getStrainratesLongterm(shz,rcv.dip(1)*pi/180,[0,20e3],[-140e3,35e3],1);
shz.e22pl = e22_dev.*Vpl;% 1/s
shz.e23pl = -e23.*Vpl;% 1/s
%% compute stress interaction kernels
% evl contains the following as N-d matrices
% KK - fault-fault interactions [rcv.N x rcv.N]
% KL - fault-shz interactions [shz.N x rcv.N x 2]
% LK - shz-fault interactions [rcv.N x shz.N x 2]
% LL - shz-shz interactions [shz.N x shz.N x 2 x 2]

% evl_orig = computeAllStressKernelsBem(rcv,shz,boundary);
% load('kernels/evl_orig.mat','evl_orig');

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

