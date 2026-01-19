function [e22_dev, e23] = getStrainratesLongterm(shz,dip,ocean_translate_params,continental_translate_params,n)
% GETSTRAINRATESLONGTERM Compute long-term strain rates in subduction zone shear zones
%
%   [e22_dev, e23] = getStrainratesLongterm(shz, dip, ocean_translate_params, ...
%                                                continental_translate_params, n)
%
%   Calculates deviatoric strain rate components for triangular mesh elements in 
%   a subduction zone by evaluating analytical corner flow velocity fields and 
%   computing spatial velocity gradients. The function handles both continental 
%   mantle wedge and oceanic sub-slab mantle regions with appropriate coordinate 
%   transformations.
%
% METHODOLOGY:
%   1. Determines whether each mesh element lies in continental or oceanic mantle
%      based on its position relative to the dipping slab interface
%   2. Applies region-specific coordinate translations to align with corner flow
%      coordinate system (origin at trench, with small offset to avoid singularity)
%   3. Evaluates corner flow velocity at triangle vertices using analytical solutions
%   4. Fits a planar velocity field to the three vertices to obtain velocity gradients
%   5. Computes strain rate tensor components from velocity gradient tensor
%
% INPUTS:
%   shz - Shear zone structure containing mesh geometry with fields:
%         .N       - Number of triangular elements
%         .xc      - [N x 2] element centroids [x, z]
%         .A,.B,.C - [N x 2] vertex coordinates for each triangle [x, z]
%
%   dip - Slab dip angle in radians (e.g., pi/4 for 45°)
%
%   ocean_translate_params - [1 x 2] translation vector [dx, dz] for oceanic mantle
%                            Shifts coordinates to corner flow reference frame
%                            (origin nominally at trench, offset by ~plate thickness)
%
%   continental_translate_params - [1 x 2] translation vector [dx, dz] for continental wedge
%                                  Shifts coordinates to avoid singularity at (0,0)
%
%   n - Power-law rheology exponent: 1 (Newtonian) or 3 (power-law viscous)
%
% OUTPUTS:
%   e22_dev - [N x 1] Deviatoric normal strain rate component ε̇_xx^dev
%             Computed as 0.5*(∂u/∂x - ∂v/∂z)
%
%   e23     - [N x 1] Shear strain rate component ε̇_xz  
%             Computed as 0.5*(∂u/∂z + ∂v/∂x)
%
% COORDINATE SYSTEM:
%   - x: horizontal (positive away from trench)
%   - z: vertical (negative downward into Earth)
%   - Slab dips at angle 'dip' from horizontal
%   - Continental mantle: above slab (z + x*tan(dip) ≥ 0)
%   - Oceanic mantle: below slab (z + x*tan(dip) < 0)
%
% DEPENDENCIES:
%   computeCornerflowVelocityfield - Analytical corner flow velocity solver
%
% NOTES:
%   - Translation parameters should be tuned to mesh geometry to avoid 
%     evaluation at or too near the corner flow singularity at (0,0)
%   - The -1 offset in z-translations approximates shifting by plate thickness
%   - Planar velocity assumption is valid for sufficiently small elements
%
% EXAMPLE:
%   dip = 45 * pi/180;  % 45-degree dipping slab
%   ocean_trans = [0, -20e3];     % oceanic lithosphere ~20 km thick
%   cont_trans = [-140e3, -35e3]; % offset from trench, ~35 km thick
%   [e22_dev, e23] = getStrainratesLongterm(shz, dip, ocean_trans, cont_trans, 1);
%
% AUTHORS:
%   Rishav Mallick (EOS), 2025
%
% SEE ALSO:
%   computeCornerflowVelocityfield

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