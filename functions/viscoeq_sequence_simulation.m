function [tfull,Yfull] = viscoeq_sequence_simulation(y_odefunction,Y0,rcv,shz,stress_change,Trecur)
% VISCOEQ_SEQUENCE_SIMULATION  Integrate one cycle with imposed earthquakes.
%
%   [tfull, Yfull] = viscoeq_sequence_simulation(y_odefunction, Y0, rcv, ...
%                                                shz, stress_change, Trecur)
%
%   Integrates the state from Y0 with ode45 from t = 0 to the first event,
%   applies that event's stress change to the state, integrates to the
%   next event, and so on until Trecur. An event adds dtau/Asigma to
%   log(v) on each fault patch and dsigma22, dsigma23 to the mantle
%   stresses.
%
% INPUTS
%   y_odefunction  function handle @(t, y) returning dY/dt, e.g.
%                  @(t, y) ode_viscoelastic(t, y, rcv, shz, evl)
%   Y0             [n x 1] initial state (layout as in ode_viscoelastic)
%   rcv            geometry.receiver object (uses .N, .dgf, .Asigma)
%   shz            geometry.shearZoneReceiver object (uses .N, .dgf)
%   stress_change  struct of imposed events
%                  .Nevents         number of events per cycle
%                  .Timing          [1 x Nevents] increasing event times
%                                   from the start of the cycle (s)
%                  .dtau            [N x Nevents] shear stress change (MPa)
%                  .dsigma22        [M x Nevents] change in s22 (MPa)
%                  .dsigma23        [M x Nevents] change in s23 (MPa)
%   Trecur         cycle length (s)
%
% OUTPUTS
%   tfull          [Nt x 1] times from the start of the cycle (s)
%   Yfull          [Nt x n] states at those times
%
% NOTES
%   - Uses ode45 (adaptive Runge-Kutta 4(5)) with AbsTol = RelTol = 1e-10,
%     InitialStep = 1e-6 s and MaxStep = 1e8 s.
%   - Each event time appears twice in tfull: the states just before and
%     just after the stress step.
%   - Events are applied by the local function apply_initialconditions.
%
% SEE ALSO
%   runImposedEarthquakeCycles, ode_viscoelastic, ode45
%
% AUTHORS
%   Rishav Mallick, 2023

options=odeset('Refine',1,'AbsTol',1e-10,'RelTol',1e-10,'InitialStep',1e-6,'MaxStep',1e8);

Nevents = stress_change.Nevents;
T_events = stress_change.Timing;

for count = 1:Nevents
    
    if count == 1
        if Nevents > 1 % for more than 1 earthquake
            Tseqend = T_events(count+1) - T_events(count);
        else % if it is a characteristic event
            Tseqend = Trecur - T_events(count);
        end
        
        % from 0 to event 1
        [t,Y]=ode45(y_odefunction,[0 T_events(count)],Y0,options);
        tfull = t;
        Yfull = Y;
        
        % from event 1 to event 2 OR Trecur
        Y0 = apply_initialconditions(Y,rcv,shz,stress_change,count);

        [t,Y]=ode45(y_odefunction,[0 Tseqend],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
        
    elseif count>1 && count<Nevents
        Y0 = apply_initialconditions(Y,rcv,shz,stress_change,count);
        
        [t,Y]=ode45(y_odefunction,[0 T_events(count+1)-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    else % for last event solve till Trecur (end of cycle)
        Y0 = apply_initialconditions(Y,rcv,shz,stress_change,count);

        [t,Y]=ode45(y_odefunction,[0 Trecur-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    end
end
        
end

function Y0 = apply_initialconditions(Y,rcv,shz,stress_change,count)
% function to apply initial conditions for a given state vector

Y0 = zeros(rcv.N*rcv.dgf + shz.N*shz.dgf,1);
logV_0 = Y(end,1 : rcv.dgf : rcv.dgf*rcv.N)' + stress_change.dtau(:,count)./rcv.Asigma;
s22_0 = Y(end,rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N)' + stress_change.dsigma22(:,count);
s23_0 = Y(end,rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N)' + stress_change.dsigma23(:,count);

Y0(1 : rcv.dgf : rcv.dgf*rcv.N) = logV_0;
Y0(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s22_0;
Y0(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s23_0;

end