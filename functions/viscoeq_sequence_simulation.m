function [tfull,Yfull] = viscoeq_sequence_simulation(y_odefunction,Y0,rcv,shz,stress_change,Trecur,nreps)
% VISCOEQ_SEQUENCE_SIMULATION  Integrate cycles with imposed earthquakes.
%
%   [tfull, Yfull] = viscoeq_sequence_simulation(y_odefunction, Y0, rcv, ...
%                                                shz, stress_change, Trecur)
%   [tfull, Yfull] = viscoeq_sequence_simulation(..., nreps)
%
%   Integrates the state from Y0 with ode45 through nreps cycles of length
%   Trecur. A cycle with K imposed events at times T1 < ... < TK has K+1
%   segments, [0,T1], [T1,T2], ..., [TK,Trecur]; each segment after the
%   first starts by applying the stress change of the event at its start.
%   An event adds dtau/Asigma to the fault state on each patch (log(v),
%   or v for rcv.faultLaw = 'linear'; see ode_viscoelastic) and dsigma22,
%   dsigma23 to the mantle stresses. Each cycle starts from the final
%   state of the previous one, and only the last cycle is returned.
%
% INPUTS
%   y_odefunction  function handle @(t, y) returning dY/dt, e.g.
%                  @(t, y) ode_viscoelastic(t, y, rcv, shz, evl)
%   Y0             [n x 1] initial state (layout as in ode_viscoelastic)
%   rcv            geometry.receiver object (uses .N, .dgf, .Asigma)
%   shz            geometry.shearZoneReceiver object (uses .N, .dgf)
%   stress_change  struct of imposed events (any number K >= 1)
%                  .Nevents         number of events per cycle, K
%                  .Timing          [1 x K] increasing event times from
%                                   the start of the cycle, 0 <= T < Trecur
%                                   (s); an event at 0 is applied to the
%                                   initial state
%                  .dtau            [N x K] shear stress change (MPa)
%                  .dsigma22        [M x K] change in s22 (MPa)
%                  .dsigma23        [M x K] change in s23 (MPa)
%   Trecur         cycle length (s)
%   nreps          number of cycles to integrate; only the last one is
%                  returned  [1]
%
% OUTPUTS
%   tfull          [Nt x 1] times from the start of the last cycle (s)
%   Yfull          [Nt x n] states at those times
%
% NOTES
%   - Uses ode45 (adaptive Runge-Kutta 4(5)) with AbsTol = RelTol = 1e-10,
%     InitialStep = 1e-6 s and MaxStep = 1e8 s, restarted on each segment.
%     The system is autonomous, so each segment is integrated from t = 0.
%   - Each event time appears twice in tfull: the states just before and
%     just after the stress step (an event at t = 0 appears once, after
%     the step).
%   - With nreps > 1, prints the cycle number at the start of each cycle.
%   - The cycle boundary is kept as a segment end (with no stress step),
%     so the returned cycle starts exactly at t = 0.
%
% SEE ALSO
%   runImposedEarthquakeCycles, ode_viscoelastic, ode45
%
% AUTHORS
%   Rishav Mallick, 2023

if nargin < 7
    nreps = 1;
end
Tk = stress_change.Timing(:);
if numel(Tk) ~= stress_change.Nevents
    error('viscoeq_sequence_simulation:events', ...
        'stress_change.Timing must have Nevents = %d entries.', stress_change.Nevents);
end
if isempty(Tk) || any(diff(Tk) <= 0) || Tk(1) < 0 || Tk(end) >= Trecur
    error('viscoeq_sequence_simulation:timing', ...
        'Event times must be strictly increasing, with 0 <= T < Trecur.');
end

options = odeset('Refine',1,'AbsTol',1e-10,'RelTol',1e-10,'InitialStep',1e-6,'MaxStep',1e8);

% every segment of every cycle, in order
tb = [0; Tk; Trecur];                          % segment ends within a cycle
[seg, rep] = ndgrid(1:numel(tb)-1, 1:nreps);

y = Y0(:);
tfull = []; Yfull = [];
for s = 1:numel(seg)
    j = seg(s);
    if j == 1 && nreps > 1
        disp(['Cycle number ' num2str(rep(s)) '/' num2str(nreps)])
    end
    if j > 1                                   % event j-1 starts segment j
        y = applyStressStep(y, rcv, shz, stress_change, j-1);
    end
    if tb(j+1) == tb(j)                        % event at t = 0: no segment
        continue
    end
    [t, Y] = ode45(y_odefunction, [0 tb(j+1)-tb(j)], y, options);
    y = Y(end,:)';
    if rep(s) == nreps                         % keep the last cycle only
        tfull = [tfull; t + tb(j)]; %#ok<AGROW>
        Yfull = [Yfull; Y];         %#ok<AGROW>
    end
end

end

function y = applyStressStep(y,rcv,shz,stress_change,k)
% add the stress change of event k to the state vector

iV   = 1 : rcv.dgf : rcv.dgf*rcv.N;
i22  = rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N;
i23  = rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N;

y(iV)  = y(iV)  + stress_change.dtau(:,k)./rcv.Asigma;
y(i22) = y(i22) + stress_change.dsigma22(:,k);
y(i23) = y(i23) + stress_change.dsigma23(:,k);

end
