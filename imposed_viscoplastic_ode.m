 function [Yp]= imposed_viscoplastic_ode(~,Y,flt,shz)

% fault slip velocity 
VS = flt.a>flt.b;
V_dummy=flt.Vo.*exp(Y(2:flt.dgf:flt.N*flt.dgf));

V = zeros(flt.N,1);
V(VS) = V_dummy(VS);

% Initiate state derivative
Yp=zeros(size(Y));  
%%
% s22 and s23; s33=-s22
% stress in shz zones
tau22dev = Y(flt.N*flt.dgf+1:shz.dgf:end);  
tau23 = Y(flt.N*flt.dgf+2:shz.dgf:end);
% tau33dev = -tau22dev;

% % strain-rate change!!!!! Linear 
% e22dot=tau22dev./shz.eta;
% e23=tau23./shz.eta;

% power law 
Q = sqrt(tau22dev.^2+tau23.^2);
shz.A=ones(shz.N, 1); % placeholder; delete and fix rheology in code
shz.n=1;
e22dot= shz.A.*Q.^(shz.n-1).*tau22dev;
e23dot= shz.A.*Q.^(shz.n-1).*tau23;
% e33dot=-e22dot


% Initiate state derivative
Yp=zeros(size(Y));  

%% Fault Interactions 

% slip velocity
Yp(1:flt.dgf:flt.N*flt.dgf)=V;

% think about this a little bit more. stresses are dev in mantle (so this is
% interaction between viscous shear in the mantle affecting the fault )
kv = flt.KK*(V-flt.Vpl) + shz.LK{1}*(e22dot-shz.e22pl) + ...
     +shz.LK{2}*(e23dot-shz.e23pl) + ...
     +shz.LK{3}*-(e22dot-shz.e22pl);

kv(~VS) = 0;
%  % stressing rates
% Yp(2:flt.dgf:flt.N*flt.dgf)=kv-flt.damping.*V;

% acceleration (rate of log(V/Vo)) 
Yp(2:flt.dgf:flt.N*flt.dgf)=(kv)./((flt.a-flt.b).*flt.sigmab + flt.damping.*V);


%% Shear Zone Interactions 
% LL is LLedited and stored as LLdev. 
s22dev = flt.KL{1}-0.5*(flt.KL{1}+flt.KL{3});
s33 = flt.KL{3}-0.5*(flt.KL{1}+flt.KL{3});
Yp(flt.N*flt.dgf+1:shz.dgf:end)=shz.LLdev{1,1}*(e22dot-shz.e22pl) ...
                                 +shz.LLdev{1,2}*(e23dot-shz.e23pl) ...
                                 +s22dev*(V-flt.Vpl);
Yp(flt.N*flt.dgf+2:shz.dgf:end)=shz.LLdev{2,1}*(e22dot-shz.e22pl) ...
                                 +shz.LLdev{2,2}*(e23dot-shz.e23pl) ...
                                 +s33*(V-flt.Vpl);

% 
% Yp(flt.N*flt.dgf+3:shz.dgf:end)=e22dot;
% Yp(flt.N*flt.dgf+4:shz.dgf:end)=e23dot;

end
