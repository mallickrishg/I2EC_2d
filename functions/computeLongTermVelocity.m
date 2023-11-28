function [Vx, Vz]=computeLongTermVelocity(X,Z,Vpl,Tplate,dip,shz,boundary,ocean_translate_params,continental_translate_params)
% Function computeLongTermVelocity calculates the long-term velocity field
% at the observation points
% given the shear zone, boundary data structures and the angle of suduction 
% This function uses compute_cornerflow_velocityfield to compute velocity field 
% in continental & oceanic mantle
% INPUTS:
% X,Z - Input/observation points
% Vpl - Loading velocity
% Tplate - Thickness of the plate
% dip - Angle of subduction in degrees
% shz - data structure containing shear zone mesh coordinates
% boundary - data structure containing the boundary file
% dip - slab dip in radians
% oceanic & continental translation parameters - provide as [x,z] translations
% 
% OUTPUTS:
% Vx,Vz - horizontal and vertical velocity components
% 
% Authors:
% Sharadha Sathiakumar (EOS), 2023 

oceanic_x_translate = ocean_translate_params(1); % x translation (relative to trench)
oceanic_z_translate = ocean_translate_params(2) - 1;% plate thickness ~ 20e3 (shift slightly)
continental_x_translate = continental_translate_params(1);%-140e3; % note: move to wedge (and not to 0,0)
continental_z_translate = continental_translate_params(2) - 1;% continental plate thickness ~ 35e3 (shift slightly)
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

% oceanic mantle 
Ox=X(inreg4)+ oceanic_x_translate;Oz=Z(inreg4) + oceanic_z_translate;
[Vxt,Vzt]= compute_cornerflow_velocityfield(dip*pi/180,Ox,Oz);
ind=find(inreg4==1);
    for i=1:length(ind)
    Vx(ind(i)) = Vxt(i,2)*Vpl ;
    Vz(ind(i)) = Vzt(i,2)*Vpl;
    end

% continental mantle    
Cx=X(inreg5) + continental_x_translate;Cz=Z(inreg5) + continental_z_translate;
[Vxt,Vzt]= compute_cornerflow_velocityfield(dip*pi/180,Cx,Cz);
ind=find(inreg5==1);

    for i=1:length(ind)
    Vx(ind(i)) = Vxt(i,1)*Vpl;
    Vz(ind(i)) = Vzt(i,1)*Vpl;
    end

end 