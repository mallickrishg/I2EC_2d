function [Ktau,Ksigma] = computetractionkernels(src,rcv)
% traction kernel computation for a given source and receiver object pair
% INPUTS
% src,rcv - objects or data structures containing fault mesh 
%           (end points and center nodes)
% 
% OUTPUTS 
% 2 x [N x N] matrices containing fault-centric traction kernels
% these kernels are for half-space only
% Ktau - traction kernel in shear direction
% Ksigma - traction kernel in fault-normal direction
% 
% Rishav Mallick, EOS, 2020

Ktau = zeros(rcv.N,src.N);
Ksigma = Ktau;
x = rcv.xc(:,1);
z = rcv.xc(:,2);

for i = 1:src.N
    m = [src.x(i,2) src.x(i,1) src.W(i) src.dip(i) 1];
    [Sxx,Sxz,Szz] = EdgeStress(m,x,z,src.earthModel.nu,src.earthModel.G);
    
    t=[Sxx.*rcv.nv(:,1)+Sxz.*rcv.nv(:,2), ...
        Sxz.*rcv.nv(:,1)+Szz.*rcv.nv(:,2)];
    
    Ktau(:,i) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);
    Ksigma(:,i) = rcv.nv(:,1).*t(:,1) + rcv.nv(:,2).*t(:,2);
end



end