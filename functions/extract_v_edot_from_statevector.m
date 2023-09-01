function [V,e22dot,e23dot] = extract_v_edot_from_statevector(Y,rcv,shz)
% function to extract fault slip rate (v) and bulk viscous strain rates
% (e22dot, e23dot) from the state vector. 
% State vector (Y) is the outputs from an ODE solve.

% megathrust fault slip rates
V = exp(Y(:,1 : rcv.dgf : rcv.dgf*rcv.N));
V(:,rcv.pinnedPosition) = 0;% locked section of fault has v = 0

% strain rates
e22dot = Y(:,rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);
e23dot = Y(:,rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N);
end