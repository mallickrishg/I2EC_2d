function [Yp] = odeViscoelastic(~,Y,rcv,shz,evl)
% ODEVISCOELASTIC  Right-hand side of the coupled fault-mantle ODEs.
%
%   Yp = odeViscoelastic(t, Y, rcv, shz, evl)
%
%   Returns dY/dt for the state Y = [fault state on each patch; s22 and
%   s23 on each mantle cell, interleaved; Kelvin strains, if any]. The
%   mantle follows a power law (the Maxwell, steady-state element),
%       edotM_ij = alpha * tau^(n-1) * s_ij,   tau = sqrt(s22^2 + s23^2),
%   in series with the p Kelvin elements of shz.etaK, shz.muK (transient
%   creep: Burgers p = 1, Andrade a Prony series; none by default),
%       s_ij = muK * epsK_ij + etaK * depsK_ij/dt,
%   and the strain rate edot_ij that loads the kernels is the sum of all
%   elements'. Each Kelvin element is stored by its spring stress
%   sK_ij = muK*epsK_ij (MPa, the scale of the mantle stresses, so one
%   ode tolerance suits both), with
%       depsK_ij/dt = (s_ij - sK_ij)/etaK,   dsK_ij/dt = muK*depsK_ij/dt,
%   after the mantle stresses, per Kelvin column k = 1..p and component:
%   [sK22(:,1); sK23(:,1); sK22(:,2); ...] (2*M*p entries; zero where
%   etaK = Inf).
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
%   Y       [n x 1] state, n = N*rcv.dgf + M*shz.dgf + 2*M*p
%   rcv     geometry.receiver object (uses .N, .dgf, .Vpl, .Asigma,
%           .pinnedPosition, .faultLaw)
%   shz     geometry.shearZoneReceiver object (uses .N, .dgf, .alpha, .n,
%           .e22pl, .e23pl, and optionally .etaK, .muK: see
%           extractKelvinElements)
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
%     ODE as the semi-analytical solvers (solveAnalyticIvp) solve, with or
%     without Kelvin elements (assembleLinearSystem). The
%     slip rate may then become negative (back-slip), which log(v) could
%     not represent.
%   - Stability: the function integrates whatever kernels it is given.
%     The default kernels of computeAllStressKernelsBem (cell-averaged,
%     certified) have no growing modes for any viscosity or power-law
%     exponent. Centroid kernels need 'kernelmodify' 1 or 2.
%
% SEE ALSO
%   runImposedEarthquakeCycles, integrateImposedCycles,
%   extractRatesFromOdeState
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

% Kelvin elements (transient creep) add their strain rates
[etaK, muK, p] = extractKelvinElements(shz);
if p > 0
    iK = rcv.dgf*rcv.N + shz.dgf*shz.N + (1:2*shz.N*p);
    sK = reshape(Y(iK), shz.N, 2, p);                  % spring stresses
    eK22dot = (s22 - reshape(sK(:,1,:), shz.N, p)) ./ etaK;
    eK23dot = (s23 - reshape(sK(:,2,:), shz.N, p)) ./ etaK;
    Yp(iK) = reshape([reshape(muK.*eK22dot, shz.N, 1, p), ...
                      reshape(muK.*eK23dot, shz.N, 1, p)], [], 1);
    e22dot = e22dot + sum(eK22dot, 2);
    e23dot = e23dot + sum(eK23dot, 2);
end

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
