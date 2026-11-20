function [Yp] = ode_viscoelastic(~,Y,rcv,shz,evl)
% ODE_VISCOELASTIC  Right-hand side of the coupled fault-mantle ODEs.
%
%   Yp = ode_viscoelastic(t, Y, rcv, shz, evl)
%
%   Returns dY/dt for the state Y = [fault state on each patch; s22 and
%   s23 on each mantle cell, interleaved]. The mantle follows a power law,
%       edot_ij = alpha * tau^(n-1) * s_ij,   tau = sqrt(s22^2 + s23^2).
%   The fault law is set by rcv.faultLaw:
%     'ratestrengthening' (default): rate-dependent friction integrated in
%         log(v), fault state log(v),  d log(v)/dt = (dtau/dt) / Asigma
%     'linear': linear viscous, tau = Asigma * v (Asigma a viscosity per
%         unit length, as in the semi-analytical solvers), fault state v,
%         dv/dt = (dtau/dt) / Asigma
%   and the stressing rates come from the kernels acting on the departures
%   from the long-term rates (j = 1 for s22, 2 for s23):
%       dtau/dt = KK*(v-Vpl) + LK1*(e22dot-e22pl) + LK2*(e23dot-e23pl)
%       dsj/dt  = KLj*(v-Vpl) + LLj1*(e22dot-e22pl) + LLj2*(e23dot-e23pl)
%
% INPUTS
%   t       time (s); unused, the system is autonomous
%   Y       [n x 1] state, n = N*rcv.dgf + M*shz.dgf
%   rcv     geometry.receiver object (uses .N, .dgf, .Vpl, .Asigma,
%           .pinnedPosition, .faultLaw)
%   shz     geometry.shearZoneReceiver object (uses .N, .dgf, .alpha, .n,
%           .e22pl, .e23pl)
%   evl     struct of stress kernels (.KK, .KL, .LK, .LL) from
%           computeAllStressKernelsBem
%
% OUTPUTS
%   Yp      [n x 1] time derivative of the state
%
% NOTES
%   - Locked patches (rcv.pinnedPosition) have v = 0 in the stress
%     calculation; their fault-state entries still evolve but are not
%     used.
%   - With the linear law and n = 1 the system is linear and is the same
%     ODE as the semi-analytical solvers (analyticSolveIVP) solve. The
%     slip rate may then become negative (back-slip), which log(v) could
%     not represent.
%   - Stability: the function integrates whatever kernels it is given.
%     The default kernels of computeAllStressKernelsBem (cell-averaged,
%     certified) have no growing modes for any viscosity or power-law
%     exponent. Centroid kernels need 'kernelmodify' 1 or 2.
%
% SEE ALSO
%   runImposedEarthquakeCycles, viscoeq_sequence_simulation,
%   extract_v_edot_from_statevector
%
% AUTHORS
%   Rishav Mallick, 2023

% Initialize state derivative
Yp = zeros(size(Y));  

% fault slip velocity
linearFault = isLinearFaultLaw(rcv);
if linearFault
    V = Y(1 : rcv.dgf : rcv.dgf*rcv.N);
else
    V = exp(Y(1 : rcv.dgf : rcv.dgf*rcv.N));
end
VW = rcv.pinnedPosition;
% locked section of fault has v = 0
V(VW) = 0;

% stress in shz zones
s22 = Y(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);
s23 = Y(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);

% power law 
sigma = sqrt(s22.^2 + s23.^2);

% strain rates
e22dot = shz.alpha.*(sigma.^(shz.n-1)).*s22;
e23dot = shz.alpha.*(sigma.^(shz.n-1)).*s23;

%% Fault Interactions
dtaudt = evl.KK*(V-rcv.Vpl) + evl.LK(:,:,1)*(e22dot - shz.e22pl) + evl.LK(:,:,2)*(e23dot - shz.e23pl);
% rate-strengthening: zeta = log(v), dzeta/dt = (dtau/dt)/Asigma
% linear viscous:     tau = Asigma*v,  dv/dt  = (dtau/dt)/Asigma
% (the same expression for the derivative of the fault state)
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = dtaudt./rcv.Asigma;

%% Shear Zone Interactions 

dsigma22dt = evl.KL(:,:,1)*(V-rcv.Vpl) + evl.LL(:,:,1,1)*(e22dot - shz.e22pl) + evl.LL(:,:,1,2)*(e23dot - shz.e23pl);
dsigma23dt = evl.KL(:,:,2)*(V-rcv.Vpl) + evl.LL(:,:,2,1)*(e22dot - shz.e22pl) + evl.LL(:,:,2,2)*(e23dot - shz.e23pl);

% store stressing rate in state derivative
Yp(rcv.N*rcv.dgf+1:shz.dgf:rcv.dgf*rcv.N+shz.dgf*shz.N) = dsigma22dt;
Yp(rcv.N*rcv.dgf+2:shz.dgf:rcv.dgf*rcv.N+shz.dgf*shz.N) = dsigma23dt;


end
