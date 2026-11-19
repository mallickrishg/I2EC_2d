function [Ktau,Ksigma] = computetractionkernels(src,rcv)
% COMPUTETRACTIONKERNELS  Half-space traction kernels due to fault slip (edge dislocations).
%
%   [Ktau, Ksigma] = computetractionkernels(src, rcv)
%
%   For unit slip on each source patch, modelled as an edge dislocation
%   with EdgeStress, computes the shear and normal traction at the receiver
%   patch centres in a homogeneous half-space. There is no boundary
%   correction.
%
% INPUTS
%   src     geometry.receiver object: source patches (Ns patches)
%           .x    [Ns x 2] upper end [x, z] of each patch (m)
%           .W    [Ns x 1] patch widths (m)
%           .dip  [Ns x 1] dips (degrees)
%           .earthModel  elastic parameters .G (MPa) and .nu
%   rcv     geometry.receiver object: receiver patches (Nr patches)
%           .xc   [Nr x 2] patch centres (m)
%           .dv   [Nr x 2] unit vectors in the dip direction
%           .nv   [Nr x 2] unit normal vectors
%
% OUTPUTS
%   Ktau    [Nr x Ns] shear traction in the receiver dip direction (MPa/m)
%   Ksigma  [Nr x Ns] normal traction (MPa/m)
%
% NOTES
%   - Not used by the I2EC solvers, which build fault kernels with
%     geometry.computeFaultTractionKernels (line dislocations, LDstressHS)
%     through computeFaultTractionKernelsBem.
%   - The +geometry package has a copy that calls geometry.EdgeStress.
%
% SEE ALSO
%   EdgeStress, computeFaultTractionKernelsBem
%
% AUTHORS
%   Rishav Mallick (EOS), 2020

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