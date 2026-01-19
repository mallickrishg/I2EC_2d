function [e22_dev, e23] = getStrainratesLongterm(shz,dip,ocean_translate_params,continental_translate_params)
% Function edotlongtermsol calculates the long-term strain rate 
% given the shear zone data structure and the angle of suduction 
% This function uses compute_cornerflow_velocityfield to compute velocity field 
% in continental & oceanic mantle
% INPUTS:
% shz - data structure containing shear zone mesh coordinates
% dip - slab dip in radians
% oceanic & continental translation parameters - provide as [x,z] translations
% 
% OUTPUTS:
% e22,e23 - strain rate in deviatoric xx and xz components
% 
% Authors:
% Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

% rescale parameters: move the shz the x and z coordinates because corner flow solutions are at (0,0)
% Note: Don't move to 0 because solutions are not well defined at 0,0. These
% parameters have to be adjusted based on how the mesh is created. 

oceanic_x_translate = ocean_translate_params(1); % x translation (relative to trench)
oceanic_z_translate = ocean_translate_params(2) - 1;% plate thickness ~ 20e3 (shift slightly)

continental_x_translate = continental_translate_params(1);%-140e3; % note: move to wedge (and not to 0,0)
continental_z_translate = continental_translate_params(2) - 1;% continental plate thickness ~ 35e3 (shift slightly)

% to calculate velocity gradient using finite difference method 
epsilonx=1e-6;
epsilonz=1e-6;

% long-term strain rates
e22_dev = zeros(shz.N,1);
e23 = zeros(shz.N,1);

for i=1:shz.N

% shear zone center
Cx=shz.xc(i,1);Cz=shz.xc(i,2);

% check whether the vertices are in the continental or oceanic mantle 
kk=Cz+Cx.*tan(dip);
        if kk>=0 
        % continental mantle 
        % Triangle co-ordinates are shifted 
        Cx=Cx + continental_x_translate;Cz=Cz + continental_z_translate;
        [Vx,Vz]     = compute_cornerflow_velocityfield(dip,Cx,Cz);
        [Vxdx,Vzdx] = compute_cornerflow_velocityfield(dip,Cx+epsilonx,Cz);
        [Vxdz,Vzdz] = compute_cornerflow_velocityfield(dip,Cx,Cz+epsilonz);
        
        e22_dev(i) = 0.5*((Vxdx(1,1)-Vx(1,1))/(epsilonx)-(Vzdz(1,1)-Vz(1,1))/(epsilonz));        
        e23(i) = 0.5*((Vxdz(1,1)-Vx(1,1))/(epsilonz)+(Vzdx(1,1)-Vz(1,1))/(epsilonx));
        % e22_dev(i) = 0.5*((Vxdx(1,1)-Vx(1,1))/(epsilonx) + (Vzdz(1,1)-Vz(1,1))/(epsilonz));
        % e23(i) = 0.5*((Vxdz(1,1)-Vx(1,1))/(epsilonz)-(Vzdx(1,1)-Vz(1,1))/(epsilonx));

        else  

        % Triangle co-ordinates are shifted 
        Cx=Cx + oceanic_x_translate;Cz=Cz + oceanic_z_translate;
        [Vx,Vz]     = compute_cornerflow_velocityfield(dip,Cx,Cz);
        [Vxdx,Vzdx] = compute_cornerflow_velocityfield(dip,Cx+epsilonx,Cz);
        [Vxdz,Vzdz] = compute_cornerflow_velocityfield(dip,Cx,Cz+epsilonz);
        
        e22_dev(i) = 0.5*((Vxdx(1,2)-Vx(1,2))/(epsilonx)-(Vzdz(1,2)-Vz(1,2))/(epsilonz));
        e23(i) = 0.5*((Vxdz(1,2)-Vx(1,2))/(epsilonz)+(Vzdx(1,2)-Vz(1,2))/(epsilonx));
        % e22_dev(i) = 0.5*((Vxdx(1,2)-Vx(1,2))/(epsilonx) - (Vzdz(1,2)-Vz(1,2))/(epsilonz));
        % e23(i) = 0.5*((Vxdz(1,2)-Vx(1,2))/(epsilonz) + (Vzdx(1,2)-Vz(1,2))/(epsilonx));

        end 


 
end 