function [Ktau,Ksigma] = computetractionkernels(src,rcv)
% COMPUTETRACTIONKERNELS  Fault tractions from slip with edge dislocations.
%
%   [Ktau, Ksigma] = geometry.computetractionkernels(src, rcv)
%
%   Shear and normal traction at the centres of rcv due to unit slip on
%   each patch of src, built from pairs of edge dislocations
%   (geometry.EdgeStress). An older alternative to
%   geometry.computeFaultTractionKernels; not called by the I2EC drivers.
%
% INPUTS
%   src     geometry.receiver object: source patches
%           .N           [1 x 1] number of patches
%           .x           [N x 2] upper ends of the patches (m, z negative)
%           .W           [N x 1] patch widths (m)
%           .dip         [N x 1] patch dips (degrees)
%           .earthModel  geometry.LDhs object (uses .G, .nu)
%   rcv     geometry.receiver object: receiver patches
%           .N           [1 x 1] number of patches
%           .xc          [Nr x 2] patch centres (m, z negative)
%           .dv          [Nr x 2] unit dip vectors
%           .nv          [Nr x 2] unit normal vectors
%
% OUTPUTS
%   Ktau    [Nr x N] shear traction per unit slip (MPa/m)
%   Ksigma  [Nr x N] normal traction per unit slip (MPa/m)
%
% SEE ALSO
%   geometry.EdgeStress, geometry.computeFaultTractionKernels
%
% AUTHORS
%   Rishav Mallick, 2020

Ktau = zeros(rcv.N,src.N);
Ksigma = Ktau;
x = rcv.xc(:,1);
z = rcv.xc(:,2);

for i = 1:src.N
    m = [src.x(i,2) src.x(i,1) src.W(i) src.dip(i) 1];
    [Sxx,Sxz,Szz] = geometry.EdgeStress(m,x,z,src.earthModel.nu,src.earthModel.G);
    
    t=[Sxx.*rcv.nv(:,1)+Sxz.*rcv.nv(:,2), ...
        Sxz.*rcv.nv(:,1)+Szz.*rcv.nv(:,2)];
    
    Ktau(:,i) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);
    Ksigma(:,i) = rcv.nv(:,1).*t(:,1) + rcv.nv(:,2).*t(:,2);
end



end