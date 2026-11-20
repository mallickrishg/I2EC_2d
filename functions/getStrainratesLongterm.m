function [e22_dev, e23] = getStrainratesLongterm(shz,dip,ocean_translate_params,continental_translate_params,n)
% GETSTRAINRATESLONGTERM  Compute long-term mantle strain rates (corner flow).
%
%   [e22_dev, e23] = getStrainratesLongterm(shz, dip, ...
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
%     The two translations move each region into the corner-flow frame and
%     keep the evaluation points away from that singularity; 1 m is also
%     subtracted from both z shifts. Tune them to the mesh.
%   - Wedge cells use column 1 (continental arc) of
%     computeCornerflowVelocityfield, oceanic cells column 2.
%   - The planar fit assumes the velocity varies little across a cell, so
%     it is accurate for cells that are small compared with their distance
%     from the corner.
%   - Rates are for a unit plate speed: the driver scripts multiply by Vpl
%     and use -e23 as the long-term shear strain rate (shz.e23pl).
%
% EXAMPLE
%   dip = 45*pi/180;
%   ocean_trans = [0, -20e3];       % oceanic lithosphere ~20 km thick
%   cont_trans = [-140e3, -35e3];   % offset from trench, ~35 km thick
%   [e22_dev, e23] = getStrainratesLongterm(shz, dip, ocean_trans, ...
%                                           cont_trans, 1);
%
% SEE ALSO
%   computeCornerflowVelocityfield, computeLongTermVelocity
%
% AUTHORS
%   Rishav Mallick (EOS), 2025

if nargin < 5
    n = 1;
end

% rescale parameters: move the shz the x and z coordinates because corner flow solutions are at (0,0)
% Note: Don't move to 0 because solutions are not well defined at 0,0. These
% parameters have to be adjusted based on how the mesh is created. 

oceanic_x_translate = ocean_translate_params(1); % x translation (relative to trench)
oceanic_z_translate = ocean_translate_params(2) - 1;% plate thickness ~ 20e3 (shift slightly)

continental_x_translate = continental_translate_params(1);%-140e3; % note: move to wedge (and not to 0,0)
continental_z_translate = continental_translate_params(2) - 1;% continental plate thickness ~ 35e3 (shift slightly)

% long-term strain rates
e22_dev = zeros(shz.N,1);
e23 = zeros(shz.N,1);

% loop through all triangles
for i = 1:shz.N

    % check whether the vertices are in the continental or oceanic mantle
    kk = shz.xc(i,2) + shz.xc(i,1).*tan(dip);

    if kk >= 0
        % continental mantle
        % Triangle co-ordinates are shifted
        Cx = [shz.A(i,1),shz.B(i,1),shz.C(i,1),shz.xc(i,1)]' + continental_x_translate;
        Cz = [shz.A(i,2),shz.B(i,2),shz.C(i,2),shz.xc(i,2)]' + continental_z_translate;
        % Velocity vectors
        [Vx,Vz]     = computeCornerflowVelocityfield(dip,Cx,Cz,n);
        Vx = Vx(:,1);
        Vz = Vz(:,1);

    else
        % oceanic mantle
        % Triangle co-ordinates are shifted
        Cx = [shz.A(i,1),shz.B(i,1),shz.C(i,1),shz.xc(i,1)]' + oceanic_x_translate;
        Cz = [shz.A(i,2),shz.B(i,2),shz.C(i,2),shz.xc(i,2)]' + oceanic_z_translate;
        % Velocity vectors
        [Vx,Vz]     = computeCornerflowVelocityfield(dip,Cx,Cz,n);
        Vx = Vx(:,2);
        Vz = Vz(:,2);
    end    

    % Build the system
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