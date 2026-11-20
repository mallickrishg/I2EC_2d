function [Vx, Vz]=computeLongTermVelocity(X,Z,Vpl,Tplate,dip,shz,boundary,ocean_translate_params,continental_translate_params)
% COMPUTELONGTERMVELOCITY  Compute the long-term velocity at observation points.
%
%   [Vx, Vz] = computeLongTermVelocity(X, Z, Vpl, Tplate, dip, shz, ...
%       boundary, ocean_translate_params, continental_translate_params)
%
%   Assigns the long-term (geological) velocity to each point by region:
%   the incoming plate moves at Vpl, the overriding plate is fixed, the
%   slab moves down-dip at Vpl, and the oceanic mantle and mantle wedge
%   follow Newtonian corner flow (cornerFlowVelocity, n = 1) scaled
%   by Vpl. Points in no region get zero velocity.
%
% INPUTS
%   X, Z                          [Nobs x 1] observation points (m), Z negative
%                                 below the surface
%   Vpl                           plate speed (m/s)
%   Tplate                        plate (lithosphere) thickness (m)
%   dip                           slab dip (degrees)
%   shz                           geometry.shearZoneReceiver object: mantle mesh
%                                 .vert     [Nv x 2] vertices [x, z] (m)
%   boundary                      geometry.receiver object: model boundary
%                                 .x        [Nb x 2] patch upper ends (m)
%   ocean_translate_params        [1 x 2] shift [dx, dz] (m) of oceanic-mantle
%                                 points into the corner-flow frame
%   continental_translate_params  [1 x 2] shift [dx, dz] (m) of mantle-wedge
%                                 points into the corner-flow frame
%
% OUTPUTS
%   Vx      [Nobs x 1] horizontal velocity (m/s)
%   Vz      [Nobs x 1] vertical velocity (m/s)
%
% NOTES
%   - Regions are polygons built from the boundary and mesh extents
%     (inpolygon): incoming plate (x < 0, above -Tplate), overriding plate,
%     slab of thickness Tplate, oceanic mantle and mantle wedge. Where they
%     overlap, the region assigned last wins (oceanic, then wedge mantle).
%   - The translations move each region into the corner-flow frame
%     (cornerFlowVelocity, which also subtracts 1 m from the z shift to
%     keep points off the corner), as in getStrainratesLongterm.
%   - dip is in degrees here (getStrainratesLongterm takes radians).
%
% SEE ALSO
%   cornerFlowVelocity, getStrainratesLongterm, inpolygon
%
% AUTHORS
%   Sharadha Sathiakumar (EOS), 2023

om=shz.vert(:,1) < -shz.vert(:,2)/tand(dip);

pt1=min(shz.vert(~om,2))+((max(shz.vert(:,1))-max(boundary.x(:,1)))/cosd(dip))*sind(dip);

 
reg1x=[min(boundary.x(:,1)) 0 0 min(boundary.x(:,1))];
reg1y=[0 0 -Tplate -Tplate];

reg2x=[0 max(boundary.x(:,1)) max(boundary.x(:,1))  min(shz.vert(~om,1)) min(shz.vert(~om,1))...
 0];
reg2y=[0 0 max(shz.vert(~om,2)/1e3) max(shz.vert(~om,2)/1e3) min(shz.vert((find(shz.vert(:,1)==min(shz.vert(~om,1)))),2))/1e3 0]*1e3;

reg3x=[0 max(boundary.x(:,1)) max(boundary.x(:,1)) 0 0];
reg3y=[0 pt1 pt1-Tplate -Tplate 0];

reg4x=[min(boundary.x(:,1)) 0 max(boundary.x(:,1)) max(boundary.x(:,1)) min(boundary.x(:,1))];
reg4y=[-Tplate -Tplate pt1-Tplate min(boundary.x(:,2)) min(boundary.x(:,2))];


reg5x=[min(shz.vert(~om,1))/1e3 min(shz.vert(~om,1))/1e3 max(boundary. x(:,1))/1e3 max(boundary.x(:,1))/1e3 min(shz.vert(~om,1))/1e3]*1e3;
reg5y=[min(shz.vert((find(shz.vert(:,1)==min(shz.vert(~om,1)))),2)) max(shz.vert(~om,2)) max(shz.vert(~om,2)) pt1  min(shz.vert((find(shz.vert(:,1)==min(shz.vert(~om,1)))),2))];

[inreg1,~]=inpolygon(X,Z,reg1x,reg1y);
[inreg2,~]=inpolygon(X,Z,reg2x,reg2y);
[inreg3,~]=inpolygon(X,Z,reg3x,reg3y);
[inreg4,~]=inpolygon(X,Z,reg4x,reg4y);
[inreg5,~]=inpolygon(X,Z,reg5x,reg5y);

Vx= zeros(size(X,1),1);
Vz= zeros(size(X,1),1);

Vx(inreg1) = Vpl;
Vz(inreg1) = 0;

Vx(inreg2) = 0;
Vz(inreg2) = 0;

Vx(inreg3) = Vpl*cosd(dip);
Vz(inreg3) = -Vpl*sind(dip);

% oceanic mantle and mantle wedge: Newtonian corner flow
[Vx(inreg4), Vz(inreg4)] = cornerFlowVelocity(dip*pi/180, X(inreg4), Z(inreg4), 'ocean', ocean_translate_params, 1);
[Vx(inreg5), Vz(inreg5)] = cornerFlowVelocity(dip*pi/180, X(inreg5), Z(inreg5), 'wedge', continental_translate_params, 1);
Vx(inreg4 | inreg5) = Vpl*Vx(inreg4 | inreg5);
Vz(inreg4 | inreg5) = Vpl*Vz(inreg4 | inreg5);

end
