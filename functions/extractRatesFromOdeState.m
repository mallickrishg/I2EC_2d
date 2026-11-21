function [V,e22dot,e23dot] = extractRatesFromOdeState(Y,rcv,shz)
% EXTRACTRATESFROMODESTATE  Slip and strain rates from ODE states.
%
%   [V, e22dot, e23dot] = extractRatesFromOdeState(Y, rcv, shz)
%
%   Converts states of the numerical integration (see odeViscoelastic) to
%   fault slip rate (v = exp(state), or v = state for rcv.faultLaw =
%   'linear'; v = 0 on locked patches), and to
%   mantle strain rates from the stresses through the power law
%       edot_ij = alpha * tau^(n-1) * s_ij,   tau = sqrt(s22^2 + s23^2)
%
% INPUTS
%   Y       [Nt x n] states, one time per row (as returned by ode45)
%   rcv     geometry.receiver object (uses .N, .dgf, .pinnedPosition,
%           .faultLaw)
%   shz     geometry.shearZoneReceiver object (uses .N, .dgf, .alpha, .n)
%
% OUTPUTS
%   V       [Nt x N] fault slip rate (m/s); zero on locked patches
%   e22dot  [Nt x M] deviatoric strain rate e22 (1/s)
%   e23dot  [Nt x M] shear strain rate e23 (1/s)
%
% SEE ALSO
%   runImposedEarthquakeCycles, odeViscoelastic
%
% AUTHORS
%   Rishav Mallick, 2023

% megathrust fault slip rates
if isLinearFaultLaw(rcv)
    V = Y(:,1 : rcv.dgf : rcv.dgf*rcv.N);
else
    V = exp(Y(:,1 : rcv.dgf : rcv.dgf*rcv.N));
end
V(:,logical(rcv.pinnedPosition)) = 0;% locked section of fault has v = 0

% stress
s22 = Y(:,rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);
s23 = Y(:,rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);

% convert from stress to strain rate
sigma = sqrt(s22.^2 + s23.^2);

A_matrix = repmat(shz.alpha',length(V(:,1)),1);
n_matrix = repmat(shz.n',length(V(:,1)),1);
e22dot = A_matrix.*(sigma.^(n_matrix-1)).*s22;
e23dot = A_matrix.*(sigma.^(n_matrix-1)).*s23;

end