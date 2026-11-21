function [e22_dev, e23] = computeLongTermStrainRates(shz,dip,ocean_translate_params,continental_translate_params,n)
% COMPUTELONGTERMSTRAINRATES  Compute long-term mantle strain rates (corner flow).
%
%   [e22_dev, e23] = computeLongTermStrainRates(shz, dip, ...
%       ocean_translate_params, continental_translate_params, n)
%
%   Evaluates the analytical corner-flow velocity field (unit plate speed)
%   at the three vertices and the centroid of every shear-zone triangle,
%   fits a planar velocity field to these four points by least squares,
%   and returns the deviatoric strain rates of that planar field. Cells
%   above the slab use the mantle-wedge (continental arc) solution, cells
%   below it the oceanic-mantle solution. The result is the long-term
%   (background) viscous strain rate that loads the earthquake-cycle model.
%
% INPUTS
%   shz                           geometry.shearZoneReceiver object:
%                                 mantle cells
%                                 .N        number of triangles M
%                                 .xc       [M x 2] centroids [x, z] (m)
%                                 .A, .B, .C  [M x 2] vertices [x, z] (m)
%   dip                           slab dip (radians)
%   ocean_translate_params        [1 x 2] shift [dx, dz] (m) added to the
%                                 coordinates of oceanic-mantle cells
%   continental_translate_params  [1 x 2] shift [dx, dz] (m) added to the
%                                 coordinates of mantle-wedge cells
%   n                             power-law exponent: 1 (Newtonian) or 3.
%                                 Optional, default 1
%
% OUTPUTS
%   e22_dev                       [M x 1] deviatoric normal strain rate,
%                                 0.5*(du/dx - dv/dz), per unit plate
%                                 speed (1/m)
%   e23                           [M x 1] shear strain rate,
%                                 0.5*(du/dz + dv/dx), per unit plate
%                                 speed (1/m)
%
% NOTES
%   - A cell is in the mantle wedge if z + x*tan(dip) >= 0 at its centroid
%     (x positive away from the trench, z negative below the surface), and
%     in the oceanic mantle otherwise.
%   - The corner-flow solutions are singular at their origin (the corner).
%     The two translations move each region into the corner-flow frame
%     (computeRegionCornerFlow, which also subtracts 1 m from the z shift to
%     keep points off the corner). Tune them to the mesh.
%   - The planar fit assumes the velocity varies little across a cell, so
%     it is accurate for cells that are small compared with their distance
%     from the corner.
%   - Rates are for a unit plate speed: the driver scripts multiply by Vpl
%     and use -e23 as the long-term shear strain rate (shz.e23pl).
%
% EXAMPLE
%   dip = 45*pi/180;
%   ocean_trans = [0, 20e3];        % oceanic lithosphere ~20 km thick
%   cont_trans = [-140e3, 35e3];    % offset from trench, ~35 km thick
%   [e22_dev, e23] = computeLongTermStrainRates(shz, dip, ocean_trans, ...
%                                           cont_trans, 1);
%
% SEE ALSO
%   computeRegionCornerFlow, computeCornerFlowVelocity,
%   computeLongTermVelocity
%
% AUTHORS
%   Rishav Mallick (EOS), 2025

if nargin < 5
    n = 1;
end

% the corner-flow frame of each region (computeRegionCornerFlow applies the
% translations, which have to be adjusted to how the mesh is made)
e22_dev = zeros(shz.N,1);
e23 = zeros(shz.N,1);
wedge = shz.xc(:,2) + shz.xc(:,1).*tan(dip) >= 0;

% loop through all triangles
for i = 1:shz.N
    % velocity at the three vertices and the centroid
    Cx = [shz.A(i,1),shz.B(i,1),shz.C(i,1),shz.xc(i,1)]';
    Cz = [shz.A(i,2),shz.B(i,2),shz.C(i,2),shz.xc(i,2)]';
    if wedge(i)
        [Vx, Vz] = computeRegionCornerFlow(dip, Cx, Cz, 'wedge', continental_translate_params, n);
    else
        [Vx, Vz] = computeRegionCornerFlow(dip, Cx, Cz, 'ocean', ocean_translate_params, n);
    end

    % Build the system (in model coordinates: the frame is a translation,
    % so the gradient is the same)
    A = [ones(length(Cx),1) Cx Cz zeros(length(Cx),3);
        zeros(length(Cx),3) ones(length(Cx),1) Cx Cz];   % 6x6
    b = [Vx; Vz];

    % Solve for coefficients
    coeff = A\b;
    a1 = coeff(2);  % du/dx
    a2 = coeff(3);  % du/dy
    b1 = coeff(5);  % dv/dx
    b2 = coeff(6);  % dv/dy
    % Strain rate tensor
    e22_dev(i) = 0.5*(a1 - b2);
    e23(i) = 0.5*(a2 + b1);

end

end