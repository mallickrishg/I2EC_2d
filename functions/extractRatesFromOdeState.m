function [V,e22dot,e23dot,sK,tK] = extractRatesFromOdeState(Y,rcv,shz)
% EXTRACTRATESFROMODESTATE  Slip and strain rates from ODE states.
%
%   [V, e22dot, e23dot] = extractRatesFromOdeState(Y, rcv, shz)
%   [V, e22dot, e23dot, sK] = extractRatesFromOdeState(Y, rcv, shz)
%   [V, e22dot, e23dot, sK, tK] = extractRatesFromOdeState(Y, rcv, shz)
%
%   Converts states of the numerical integration (see odeViscoelastic) to
%   fault slip rate (v = exp(state), or v = state for rcv.faultLaw =
%   'linear'; v = 0 on locked patches), and to
%   mantle strain rates from the stresses through the power law
%       edot_ij = alpha * tau^(n-1) * s_ij,   tau = sqrt(s22^2 + s23^2)
%   plus, for cells with Kelvin elements (shz.etaK, shz.muK), their rates
%   (s_ij - sK_ij)/etaK, sK = muK*epsK the spring stress, so the strain
%   rates are totals. Likewise the slip rate includes the fault Kelvin
%   elements (rcv.etaK, rcv.muK), (Asigma*state - tK)/etaK.
%
% INPUTS
%   Y       [Nt x n] states, one time per row (as returned by ode45)
%   rcv     geometry.receiver object (uses .N, .dgf, .pinnedPosition,
%           .faultLaw, .Asigma, and .etaK, .muK if set)
%   shz     geometry.shearZoneReceiver object (uses .N, .dgf, .alpha, .n,
%           and .etaK, .muK if set)
%
% OUTPUTS
%   V       [Nt x N] fault slip rate (m/s); zero on locked patches
%   e22dot  [Nt x M] deviatoric strain rate e22 (1/s)
%   e23dot  [Nt x M] shear strain rate e23 (1/s)
%   sK      [Nt x M x 2 x p] Kelvin spring stresses muK*epsK (MPa;
%           component 22, 23 in dim 3); empty without Kelvin elements
%   tK      [Nt x N x pf] fault Kelvin spring stresses muK*sK (MPa); empty
%           without fault Kelvin elements
%
% SEE ALSO
%   runImposedEarthquakeCycles, odeViscoelastic
%
% AUTHORS
%   Rishav Mallick, 2023

% megathrust fault slip rates
zeta = Y(:,1 : rcv.dgf : rcv.dgf*rcv.N);
if isLinearFaultLaw(rcv)
    V = zeta;
else
    V = exp(zeta);
end
% fault Kelvin elements (transient afterslip) add their slip rates
[etaKf, ~, pf] = extractKelvinElements(rcv);
[~, ~, p] = extractKelvinElements(shz);
tK = [];
if pf > 0
    Nt = size(Y, 1);
    iKf = rcv.dgf*rcv.N + shz.dgf*shz.N + 2*shz.N*p + (1:rcv.N*pf);
    tK = reshape(Y(:,iKf), Nt, rcv.N, pf);
    for k = 1:pf
        V = V + (rcv.Asigma(:)'.*zeta - tK(:,:,k)) ./ etaKf(:,k)';
    end
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

% Kelvin elements (transient creep) add their strain rates
[etaK, ~, p] = extractKelvinElements(shz);
sK = [];
if p > 0
    Nt = size(Y, 1); M = shz.N;
    iK = rcv.dgf*rcv.N + shz.dgf*M + (1:2*M*p);
    sK = reshape(Y(:,iK), Nt, M, 2, p);
    for k = 1:p
        e22dot = e22dot + (s22 - sK(:,:,1,k)) ./ etaK(:,k)';
        e23dot = e23dot + (s23 - sK(:,:,2,k)) ./ etaK(:,k)';
    end
end

end