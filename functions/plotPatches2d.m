function plotPatches2d(rcv,toplot)
% PLOTPATCHES2D  Plot fault or boundary patches as line segments.
%
%   plotPatches2d(rcv)
%   plotPatches2d(rcv, toplot)
%
%   Draws every patch of a geometry.receiver object as a line from its
%   upper end down-dip over its width, in km. With one input the patches
%   are thick black lines; with a value per patch they are coloured by
%   that value (cline). Calls "hold on".
%
% INPUTS
%   rcv     geometry.receiver object: fault or boundary patches
%           .N        number of patches N
%           .x        [N x 2] upper end of each patch [x, z] (m)
%           .W        [N x 1] patch widths (m)
%           .dip      [N x 1] dip (degrees)
%   toplot  [N x 1] value per patch for the colour (optional)
%
% OUTPUTS
%   none (draws into the current axes)
%
% SEE ALSO
%   plotShz2d, cline
%
% AUTHORS
%   Rishav Mallick, 2023

km2m = 1e3;

xtop = rcv.x(:,1);
xbot = rcv.x(:,1) + rcv.W.*cosd(rcv.dip);

ztop = -rcv.x(:,2);
zbot = -rcv.x(:,2) + rcv.W.*sind(rcv.dip);

if nargin==2
    for i = 1:rcv.N
        cline([xtop(i) xbot(i)]./km2m,-[ztop(i) zbot(i)]./km2m,toplot(i).*[1 1]), hold on
    end
else
    for i = 1:rcv.N
        plot([xtop(i) xbot(i)]./km2m,-[ztop(i) zbot(i)]./km2m,'k-','LineWidth',2), hold on
    end
end

end