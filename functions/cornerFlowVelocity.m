function [Vx, Vz] = cornerFlowVelocity(dip, x, z, region, shift, n)
% CORNERFLOWVELOCITY  Corner-flow velocity of one mantle region at model points.
%
%   [Vx, Vz] = cornerFlowVelocity(dip, x, z, region, shift, n)
%
%   Moves model points into the corner-flow frame of one region (mantle
%   wedge or oceanic mantle) and returns that region's corner-flow
%   velocity there, for a unit plate speed. This is the one place where
%   the translate parameters of getStrainratesLongterm and
%   computeLongTermVelocity are applied.
%
% INPUTS
%   dip     slab dip (radians)
%   x, z    [N x 1] model coordinates (m), z negative below the surface
%   region  'wedge' (mantle wedge, continental arc) or 'ocean' (oceanic
%           mantle)
%   shift   [1 x 2] translation [dx, dz] (m) of that region's points into
%           the corner-flow frame (the ocean_translate_params or
%           continental_translate_params of the callers)
%   n       power-law exponent: 1 (Newtonian) or 3. Optional, default 1
%
% OUTPUTS
%   Vx      [N x 1] horizontal velocity per unit plate speed
%   Vz      [N x 1] vertical velocity per unit plate speed
%
% NOTES
%   - The points are moved to (x + dx, z + dz - 1): the corner-flow
%     solution is singular at its corner, and the extra 1 m keeps points
%     on the plate base off it.
%   - Uses computeCornerflowVelocityfield: column 1 for the wedge, column
%     2 for the oceanic mantle.
%
% SEE ALSO
%   computeCornerflowVelocityfield, getStrainratesLongterm,
%   computeLongTermVelocity

if nargin < 6
    n = 1;
end
switch region
    case 'wedge', col = 1;
    case 'ocean', col = 2;
    otherwise
        error('cornerFlowVelocity:region', ...
            'region must be ''wedge'' or ''ocean'', not ''%s''', region);
end

[Vx, Vz] = computeCornerflowVelocityfield(dip, x(:) + shift(1), z(:) + shift(2) - 1, n);
Vx = Vx(:, col);
Vz = Vz(:, col);
end
