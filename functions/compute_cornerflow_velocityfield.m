function [Vx,Vz] = compute_cornerflow_velocityfield(dip,x,z)
% Newtonian fluid corner flow solutions subject to 2 different sets of
% boundary conditions: 
% (1) continental arc
% (2) oceanic mantle
% INPUTS:
% dip - slab dip angle in radians
% x,z - normalized distance from slab (z is negative inside the Earth)
%       provide them as [Nx1] vectors
% OUTPUTS:
% Vx,Vz - corner flow velocity field as [Nx2] matrices for each component
% Vx(:,1) - continental arc solution
% Vx(:,2) - oceanic mantle solution
% similarly for Vz
% 
% Rishav Mallick, JPL, 2023

constants_arc = ...
    [0,...
    ...
    (-1).*(2.*cos(dip)+(-2).*cos(dip).*cos(2.*dip)+(-4).*...
    dip.*sin(dip)+(-2).*sin(dip).*sin(2.*dip)).*(1+(-4).*dip.^2+ ...
    (-2).*cos(2.*dip)+cos(2.*dip).^2+sin(2.*dip).^2).^(-1),...
    ...
    (2.*cos(dip)+(-2).*cos(dip).*cos(2.*dip)+(-4).*dip.*sin(dip)+( ...
    -2).*sin(dip).*sin(2.*dip)).*(1+(-4).*dip.^2+(-2).*cos(2.* ...
    dip)+cos(2.*dip).^2+sin(2.*dip).^2).^(-1),...
    ...
    (2.*dip.*csc((1/2) ...
    .*dip+(-1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2+csc((1/2) ...
    .*dip+(-1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2.*sin(2.* ...
    dip)).^(-1).*(2.*cos(dip).*csc((1/2).*dip+(-1/4).*pi).^2.* ...
    csc((1/2).*dip+(1/4).*pi).^2+(-1).*csc((1/2).*dip+(-1/4).* ...
    pi).^2.*csc((1/2).*dip+(1/4).*pi).^2.*(2.*cos(dip)+(-2).* ...
    cos(dip).*cos(2.*dip)+(-4).*dip.*sin(dip)+(-2).*sin(dip).* ...
    sin(2.*dip)).*(1+(-4).*dip.^2+(-2).*cos(2.*dip)+cos(2.*dip) ...
    .^2+sin(2.*dip).^2).^(-1)+cos(2.*dip).*csc((1/2).*dip+(-1/4) ...
    .*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2.*(2.*cos(dip)+(-2).* ...
    cos(dip).*cos(2.*dip)+(-4).*dip.*sin(dip)+(-2).*sin(dip).* ...
    sin(2.*dip)).*(1+(-4).*dip.^2+(-2).*cos(2.*dip)+cos(2.*dip) ...
    .^2+sin(2.*dip).^2).^(-1))];

constants_ocean=...
    [pi.*(2+(-2).*cos(dip)+(-2).*cos(2.*dip)+2.*cos(dip).*cos( ...
    2.*dip)+4.*dip.*sin(dip)+(-4).*pi.*sin(dip)+2.*sin(dip).* ...
    sin(2.*dip)).*((-1)+4.*dip.^2+(-8).*dip.*pi+4.*pi.^2+2.*cos( ...
    2.*dip)+(-1).*cos(2.*dip).^2+(-1).*sin(2.*dip).^2).^(-1), ...
    ...
    (-1)+(-1).*(2+(-2).*cos(dip)+(-2).*cos(2.*dip)+2.*cos(dip).* ...
    cos(2.*dip)+4.*dip.*sin(dip)+(-4).*pi.*sin(dip)+2.*sin(dip) ...
    .*sin(2.*dip)).*((-1)+4.*dip.^2+(-8).*dip.*pi+4.*pi.^2+2.* ...
    cos(2.*dip)+(-1).*cos(2.*dip).^2+(-1).*sin(2.*dip).^2).^(-1)+...
    pi.*((-2).*dip.*csc((1/2).*dip+(-1/4).*pi).^2.*csc((1/2).* ...
    dip+(1/4).*pi).^2+2.*pi.*csc((1/2).*dip+(-1/4).*pi).^2.*csc( ...
    (1/2).*dip+(1/4).*pi).^2+(-1).*csc((1/2).*dip+(-1/4).*pi) ...
    .^2.*csc((1/2).*dip+(1/4).*pi).^2.*sin(2.*dip)).^(-1).*(2.* ...
    csc((1/2).*dip+(-1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2+ ...
    (-2).*cos(dip).*csc((1/2).*dip+(-1/4).*pi).^2.*csc((1/2).* ...
    dip+(1/4).*pi).^2+csc((1/2).*dip+(-1/4).*pi).^2.*csc((1/2).* ...
    dip+(1/4).*pi).^2.*(2+(-2).*cos(dip)+(-2).*cos(2.*dip)+2.* ...
    cos(dip).*cos(2.*dip)+4.*dip.*sin(dip)+(-4).*pi.*sin(dip)+ ...
    2.*sin(dip).*sin(2.*dip)).*((-1)+4.*dip.^2+(-8).*dip.*pi+4.* ...
    pi.^2+2.*cos(2.*dip)+(-1).*cos(2.*dip).^2+(-1).*sin(2.*dip) ...
    .^2).^(-1)+(-1).*cos(2.*dip).*csc((1/2).*dip+(-1/4).*pi) ...
    .^2.*csc((1/2).*dip+(1/4).*pi).^2.*(2+(-2).*cos(dip)+(-2).* ...
    cos(2.*dip)+2.*cos(dip).*cos(2.*dip)+4.*dip.*sin(dip)+(-4).* ...
    pi.*sin(dip)+2.*sin(dip).*sin(2.*dip)).*((-1)+4.*dip.^2+(-8) ...
    .*dip.*pi+4.*pi.^2+2.*cos(2.*dip)+(-1).*cos(2.*dip).^2+(-1) ...
    .*sin(2.*dip).^2).^(-1)),...
    ...
    (2+(-2).*cos(dip)+(-2).*cos(2.*dip) ...
    +2.*cos(dip).*cos(2.*dip)+4.*dip.*sin(dip)+(-4).*pi.*sin( ...
    dip)+2.*sin(dip).*sin(2.*dip)).*((-1)+4.*dip.^2+(-8).*dip.* ...
    pi+4.*pi.^2+2.*cos(2.*dip)+(-1).*cos(2.*dip).^2+(-1).*sin( ...
    2.*dip).^2).^(-1),...
    ...
    ((-2).*dip.*csc((1/2).*dip+(-1/4).*pi) ...
    .^2.*csc((1/2).*dip+(1/4).*pi).^2+2.*pi.*csc((1/2).*dip+( ...
    -1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2+(-1).*csc((1/2) ...
    .*dip+(-1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2.*sin(2.* ...
    dip)).^(-1).*(2.*csc((1/2).*dip+(-1/4).*pi).^2.*csc((1/2).* ...
    dip+(1/4).*pi).^2+(-2).*cos(dip).*csc((1/2).*dip+(-1/4).*pi) ...
    .^2.*csc((1/2).*dip+(1/4).*pi).^2+csc((1/2).*dip+(-1/4).*pi) ...
    .^2.*csc((1/2).*dip+(1/4).*pi).^2.*(2+(-2).*cos(dip)+(-2).* ...
    cos(2.*dip)+2.*cos(dip).*cos(2.*dip)+4.*dip.*sin(dip)+(-4).* ...
    pi.*sin(dip)+2.*sin(dip).*sin(2.*dip)).*((-1)+4.*dip.^2+(-8) ...
    .*dip.*pi+4.*pi.^2+2.*cos(2.*dip)+(-1).*cos(2.*dip).^2+(-1) ...
    .*sin(2.*dip).^2).^(-1)+(-1).*cos(2.*dip).*csc((1/2).*dip+( ...
    -1/4).*pi).^2.*csc((1/2).*dip+(1/4).*pi).^2.*(2+(-2).*cos( ...
    dip)+(-2).*cos(2.*dip)+2.*cos(dip).*cos(2.*dip)+4.*dip.*sin( ...
    dip)+(-4).*pi.*sin(dip)+2.*sin(dip).*sin(2.*dip)).*((-1)+4.* ...
    dip.^2+(-8).*dip.*pi+4.*pi.^2+2.*cos(2.*dip)+(-1).*cos(2.* ...
    dip).^2+(-1).*sin(2.*dip).^2).^(-1))];

% calculate velocity field for continental arc
A = constants_arc(1);
B = constants_arc(2);
C = constants_arc(3);
D = constants_arc(4);

vx_arc = (-1).*B + (-1).*x.*(C.*x+D.*z).*(x.^2+z.^2).^(-1) + (-1).*D.*atan2(z, x);
vz_arc = A + (-1).*z.*(C.*x+D.*z).*(x.^2+z.^2).^(-1) + C.*atan2(z, x);

% calculate velocity field for continental arc
A = constants_ocean(1);
B = constants_ocean(2);
C = constants_ocean(3);
D = constants_ocean(4);

vx_ocean = (-1).*B + (-1).*x.*(C.*x+D.*z).*(x.^2+z.^2).^(-1) + (-1).*D.*atan2(z, x);
vz_ocean = A + (-1).*z.*(C.*x+D.*z).*(x.^2+z.^2).^(-1) + C.*atan2(z, x);

% store in a single matrix and return
Vx = [vx_arc,vx_ocean];
Vz = [vz_arc,vz_ocean];

end







