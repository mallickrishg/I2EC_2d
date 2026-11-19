function plotshz2d(shz,toplot)
% PLOTSHZ2D  Plot shear-zone triangles, filled by value or as outlines.
%
%   plotshz2d(shz)
%   plotshz2d(shz, toplot)
%
%   With a value per triangle, draws the mesh as filled patches coloured by
%   that value. With one input, draws each triangle outline in red (and
%   calls "hold on"). Coordinates are plotted in km.
%
% INPUTS
%   shz     geometry.shearZoneReceiver object: mantle mesh
%           .A, .B, .C  [M x 2] triangle vertices [x, z] (m)
%           .vert     [Nv x 2] vertices (m), used for outlines
%           .tri      [M x 3] vertex indices, used for outlines
%   toplot  [M x 1] value per triangle for the colour (optional)
%
% OUTPUTS
%   none (draws into the current axes)
%
% SEE ALSO
%   plotpatch2d, patch
%
% AUTHORS
%   Rishav Mallick, 2023

km2m = 1e3;

if nargin == 2
    x = [shz.A(:,1),shz.B(:,1),shz.C(:,1)]';
    y = [shz.A(:,2),shz.B(:,2),shz.C(:,2)]';
    
    patch(x./km2m,y./km2m,toplot)
else
    for i = 1:shz.N
        plot([shz.vert(shz.tri(i,:),1);shz.vert(shz.tri(i,1),1)]./km2m,...
            [shz.vert(shz.tri(i,:),2);shz.vert(shz.tri(i,1),2)]./km2m,'r-'), hold on
    end
end
end