function [tfull,Yfull] = viscoeq_sequence_simulation(y_odefunction,Y0,rcv,shz,stress_change,Trecur)

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
        Y0 = zeros(rcv.N*rcv.dgf + shz.N*shz.dgf,1);
        logV_0 = Y(end,1 : rcv.dgf : rcv.dgf*rcv.N)' + stress_change.dtau(:,count)./rcv.Asigma;
        s22_0 = Y(end,rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N)' + stress_change.dsigma22(:,count);
        s23_0 = Y(end,rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N)' + stress_change.dsigma23(:,count);
                
        Y0(1 : rcv.dgf : rcv.dgf*rcv.N) = logV_0;
        Y0(rcv.dgf*rcv.N+1 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s22_0;
        Y0(rcv.dgf*rcv.N+2 : shz.dgf : rcv.dgf*rcv.N+shz.dgf*shz.N) = s23_0;

        [t,Y]=ode45(y_odefunction,[0 Tseqend],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
        
    elseif count>1 && count<Nevents
        Yval = exp(Y(end,:)'); 
        Y0 = log(rcv.A) + rcv.n.*(log(tau_cs(:,count) + (Yval./rcv.A).^(1./rcv.n)));
        [t,Y]=ode45(y_odefunction,[0 T_events(count+1)-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    else
        Yval = exp(Y(end,:)'); 
        Y0 = log(rcv.A) + rcv.n.*(log(tau_cs(:,count) + (Yval./rcv.A).^(1./rcv.n)));
        [t,Y]=ode45(y_odefunction,[0 Trecur-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    end
end
        
end