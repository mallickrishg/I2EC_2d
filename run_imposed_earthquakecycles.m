function [t,V,e22dot,e23dot] = run_imposed_earthquakecycles(rcv,shz,evl,stress_change,nreps)
% earthquake cycle function that solves a system of coupled ODEs to compute
% time-dependent fault slip rate and deviatoric strain rates (in shear
% zones) for a given 'rcv' and 'shz' objects.
% 
% INPUTS - 
% rcv, shz      - geometry objects
% evl           - data structure that contains all the necessary stress interaction
%                 kernels
% stress_change - data structure that constains coseismic stress changes
%                 and timings of all events that repeat cyclically
% nreps         - number of repeats of periodic cycles (needed to spin up)
% AUTHORS:
% Rishav Mallick (Caltech Seismolab), 2023

% create initial values as perturbations about steady-state values
Y0 = zeros(rcv.N*rcv.dgf + shz.N*shz.dgf,1);
logV_0 = log(rcv.Vpl);
e22dot_0 = shz.e22pl;
e23dot_0 = shz.e23pl;

Y0(1 : rcv.dgf : rcv.dgf*rcv.N) = logV_0;
Y0(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = e22dot_0;
Y0(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = e23dot_0;


% start simulation (nrep-1 spin-up cycles)
[~,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl),Y0,rcv,shz,stress_change);

for i = 2:nreps-1
    Y0 = Y(end,:)';
    [~,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl),Y0,rcv,shz,stress_change);
end
Y0 = Y(end,:)';
[t,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl),Y0,rcv,shz,stress_change);

% extract important quantities from state vector
[V,e22dot,e23dot] = extract_v_edot_from_statevector(Y,rcv,shz);

end
