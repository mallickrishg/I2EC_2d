function [V,e22dot,e23dot] = extract_v_edot_from_statevector(Y,rcv,shz)
% function to extract fault slip rate (v) and bulk viscous strain rates
% (e22dot, e23dot) from the state vector. 
% State vector (Y) is the outputs from an ODE solve.

% megathrust fault slip rates
V = exp(Y(:,1 : rcv.dgf : rcv.dgf*rcv.N));
V(:,rcv.pinnedPosition) = 0;% locked section of fault has v = 0

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