function [t,V,e22dot,e23dot] = runImposedEarthquakeCycles(rcv,shz,evl,stress_change,nreps,Trecur)
% RUNIMPOSEDEARTHQUAKECYCLES  Spin up and integrate imposed earthquake cycles.
%
%   [t, V, e22dot, e23dot] = runImposedEarthquakeCycles(rcv, shz, evl, ...
%                                         stress_change, nreps, Trecur)
%
%   Integrates the coupled fault and mantle equations (ode_viscoelastic)
%   through nreps repeats of a cycle of length Trecur, with the coseismic
%   stress changes in stress_change imposed at their event times
%   (viscoeq_sequence_simulation). Each cycle starts from the final state
%   of the previous one; the first starts from the long-term state
%   (v = rcv.Vpl and the mantle stresses that give the long-term strain
%   rates). The last cycle is returned. This is the numerical path for
%   general (power-law) rheology.
%
% INPUTS
%   rcv            geometry.receiver object: megathrust patches
%                  .N               number of patches
%                  .Vpl             [N x 1] long-term slip rate (m/s)
%                  .Asigma          [N x 1] friction parameter (a-b)*sigma
%                                   (MPa)
%                  .pinnedPosition  [N x 1] logical, true for locked
%                                   patches
%                  .dgf             degrees of freedom per patch (1)
%   shz            geometry.shearZoneReceiver object: mantle cells
%                  .N               number of cells
%                  .alpha           [M x 1] fluidity, 1/viscosity
%                                   (1/(MPa s))
%                  .n               [M x 1] power-law exponent
%                  .e22pl, .e23pl   [M x 1] long-term strain rates (1/s)
%                  .dgf             degrees of freedom per cell (2)
%   evl            struct of stress kernels (.KK, .KL, .LK, .LL) from
%                  computeAllStressKernelsBem
%   stress_change  struct of imposed events (any number K >= 1, see
%                  viscoeq_sequence_simulation)
%                  .Nevents         number of events per cycle
%                  .Timing          [1 x Nevents] event times from the
%                                   start of the cycle (s)
%                  .dtau            [N x Nevents] shear stress change (MPa)
%                  .dsigma22        [M x Nevents] change in s22 (MPa)
%                  .dsigma23        [M x Nevents] change in s23 (MPa)
%   nreps          number of cycles to run (spin-up plus the returned one)
%   Trecur         cycle (supercycle) length (s)
%
% OUTPUTS
%   t              [Nt x 1] times within the last cycle (s)
%   V              [Nt x N] fault slip rate (m/s); zero on locked patches
%   e22dot         [Nt x M] deviatoric strain rate e22 (1/s)
%   e23dot         [Nt x M] shear strain rate e23 (1/s)
%
% NOTES
%   - Stability: pass the default kernels of computeAllStressKernelsBem
%     (cell-averaged and certified: no growing modes for any viscosity or
%     power-law exponent). Centroid kernels ('cellaverage', 1) have
%     spurious growing modes and need 'kernelmodify' 1 or 2.
%   - Spin-up is by brute force: nreps must be large enough for the cycle
%     to become periodic; compare successive cycles to check.
%   - nreps = 1 runs and returns a single cycle from the long-term state,
%     with no spin-up.
%
% SEE ALSO
%   viscoeq_sequence_simulation, ode_viscoelastic,
%   extract_v_edot_from_statevector, computeAllStressKernelsBem
%
% AUTHORS
%   Rishav Mallick (Caltech Seismolab), 2023

% create initial values as perturbations about steady-state values
Y0 = zeros(rcv.N*rcv.dgf + shz.N*shz.dgf,1);
% fault state: log(v), or v for the linear viscous fault law
if isLinearFaultLaw(rcv)
    logV_0 = rcv.Vpl;
else
    logV_0 = log(rcv.Vpl);
end

% TODO - change from strain rate to stress
edot_pl = sqrt(shz.e22pl.^2 + shz.e23pl.^2);
s22_0 = (edot_pl./shz.alpha).^(1./shz.n-1).*shz.e22pl./shz.alpha;
s23_0 = (edot_pl./shz.alpha).^(1./shz.n-1).*shz.e23pl./shz.alpha;

Y0(1 : rcv.dgf : rcv.dgf*rcv.N) = logV_0;
Y0(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s22_0;
Y0(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s23_0;


% nreps-1 spin-up cycles, then the returned one
[t,Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl), ...
    Y0,rcv,shz,stress_change,Trecur,nreps);

% extract important quantities from state vector
[V,e22dot,e23dot] = extract_v_edot_from_statevector(Y,rcv,shz);

end
