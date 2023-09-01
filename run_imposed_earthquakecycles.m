function [t,V,e22dot,e33dot] = run_imposed_earthquakecycles(odefunc,rcv,shz,evl,stress_change,nreps)

[~,Y] = eqsequence_viscosims(@(t,y) ode_faultviscopower(t,y,ss,evl),ss,tau_cs,Y0,length(T_events),T_events,Teq);

for i = 2:nreps-1
    v0 = exp(Y(end,:)');v0(ss.pin) = 0;
    Y0 = log(ss.A.*((v0./ss.A).^(1./ss.n) + tau_cs(:,1)).^ss.n);
    [~,Y] = eqsequence_viscosims(@(t,y) ode_faultviscopower(t,y,ss,evl),ss,tau_cs,Y0,length(T_events),T_events,Teq);
end

v0 = exp(Y(end,:)');v0(ss.pin) = 0;
Y0 = log(ss.A.*((v0./ss.A).^(1./ss.n) + tau_cs(:,1)).^ss.n);
[tsim,Ysim] = eqsequence_viscosims(@(t,y) ode_faultviscopower(t,y,ss,evl),ss,tau_cs,Y0,length(T_events),T_events,tmodel(end));


end

function [tfull,Yfull] = eqsequence_viscosims(yp,rcv,tau_cs,Y0,Nevents,T_events,Tend)

options=odeset('Refine',1,'AbsTol',1e-10,'RelTol',1e-10,'InitialStep',1e-6,'MaxStep',1e8);

for count = 1:Nevents
    
    if count == 1
        if length(T_events)>1
            Tseqend = T_events(count+1);
        else
            Tseqend = Tend;
        end
        [t,Y]=ode45(yp,[0 Tseqend],Y0,options);
        tfull = t;
        Yfull = Y;
        
    elseif count>1 && count<Nevents
        Yval = exp(Y(end,:)'); 
        Y0 = log(rcv.A) + rcv.n.*(log(tau_cs(:,count) + (Yval./rcv.A).^(1./rcv.n)));
        [t,Y]=ode45(yp,[0 T_events(count+1)-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    else
        Yval = exp(Y(end,:)'); 
        Y0 = log(rcv.A) + rcv.n.*(log(tau_cs(:,count) + (Yval./rcv.A).^(1./rcv.n)));
        [t,Y]=ode45(yp,[0 Tend-T_events(count)],Y0,options);
        tfull = [tfull;t + T_events(count)];
        Yfull = [Yfull;Y];
    
    end
end
        
end