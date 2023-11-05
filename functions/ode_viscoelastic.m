function [Yp] = ode_viscoelastic(~,Y,rcv,shz,evl)

% Initialize state derivative
Yp = zeros(size(Y));  

% fault slip velocity 
V = exp(Y(1 : rcv.dgf : rcv.dgf*rcv.N));
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
% zeta = log(v)
% dzeta/dt = 1/v dv/dt (this is what we integrate)
dtaudt = evl.KK*(V-rcv.Vpl) + evl.LK(:,:,1)*(e22dot - shz.e22pl) + evl.LK(:,:,2)*(e23dot - shz.e23pl);
dzetadt = dtaudt./rcv.Asigma;

% d(log(v))/dt in state derivative
Yp(1:rcv.dgf:rcv.N*rcv.dgf) = dzetadt;

%% Shear Zone Interactions 

dsigma22dt = evl.KL(:,:,1)*(V-rcv.Vpl) + evl.LL(:,:,1,1)*(e22dot - shz.e22pl) + evl.LL(:,:,1,2)*(e23dot - shz.e23pl);
dsigma23dt = evl.KL(:,:,2)*(V-rcv.Vpl) + evl.LL(:,:,2,1)*(e22dot - shz.e22pl) + evl.LL(:,:,2,2)*(e23dot - shz.e23pl);

% store stressing rate in state derivative
Yp(rcv.N*rcv.dgf+1:shz.dgf:rcv.dgf*rcv.N+shz.dgf*shz.N) = dsigma22dt;
Yp(rcv.N*rcv.dgf+2:shz.dgf:rcv.dgf*rcv.N+shz.dgf*shz.N) = dsigma23dt;


end
