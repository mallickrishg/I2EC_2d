function [t,V,e22dot,e23dot] = runImposedEarthquakeCycles(rcv,shz,evl,stress_change,nreps,Trecur)
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
% Trecur        - recurrence interval for earthquake supercycle (in s)
% 
% AUTHORS:
% Rishav Mallick (Caltech Seismolab), 2023

% create initial values as perturbations about steady-state values
Y0 = zeros(rcv.N*rcv.dgf + shz.N*shz.dgf,1);
logV_0 = log(rcv.Vpl);

% TODO - change from strain rate to stress
edot_pl = sqrt(shz.e22pl.^2 + shz.e23pl.^2);
s22_0 = (edot_pl./shz.alpha).^(1./shz.n-1).*shz.e22pl./shz.alpha;
s23_0 = (edot_pl./shz.alpha).^(1./shz.n-1).*shz.e23pl./shz.alpha;

Y0(1 : rcv.dgf : rcv.dgf*rcv.N) = logV_0;
Y0(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s22_0;
Y0(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s23_0;


% start simulation (nrep-1 spin-up cycles)
[~,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl),Y0,rcv,shz,stress_change,Trecur);

for i = 2:nreps
    Y0 = Y(end,:)';
    [t,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl),Y0,rcv,shz,stress_change,Trecur);
end

% extract important quantities from state vector
[V,e22dot,e23dot] = extract_v_edot_from_statevector(Y,rcv,shz);

end
