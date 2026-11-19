function [Vx, Vz] = computeCornerflowVelocityfield(dip, x, z, n)
% COMPUTECORNERFLOWVELOCITYFIELD  Evaluate analytical subduction corner flow.
%
%   [Vx, Vz] = computeCornerflowVelocityfield(dip, x, z, n)
%
%   Returns the corner-flow velocity for a unit plate speed in two regions
%   at once: the mantle wedge (continental arc) between the fixed upper
%   plate and the slab, and the oceanic mantle below the incoming plate and
%   the slab. For n = 1 (Newtonian) the closed-form stream-function
%   solution is used; for n = 3 the coefficients of the power-law solution
%   are found from the boundary conditions.
%
% INPUTS
%   dip     slab dip (radians)
%   x       [N x 1] horizontal coordinates in the corner-flow frame
%   z       [N x 1] vertical coordinates, negative inside the Earth
%   n       power-law exponent: 1 (Newtonian) or 3
%
% OUTPUTS
%   Vx      [N x 2] horizontal velocity per unit plate speed;
%           column 1 mantle wedge (arc), column 2 oceanic mantle
%   Vz      [N x 2] vertical velocity per unit plate speed;
%           column 1 mantle wedge (arc), column 2 oceanic mantle
%
% NOTES
%   - Both columns are evaluated at every point; the caller picks the column
%     for the region each point lies in.
%   - The solution is singular at the corner (origin). Callers translate the
%     coordinates of each region (ocean/continental translate parameters in
%     getStrainratesLongterm) to keep points away from it; r = 0 is set to
%     r^2 = Inf in the Newtonian branch.
%   - Wedge: fixed upper surface, slab moving down-dip. Oceanic mantle:
%     incoming plate moving at unit speed, slab moving down-dip.
%   - For n = 3 the four coefficients of each region come from the
%     boundary conditions: a scalar equation for the phase is solved
%     with fzero and the other three follow linearly. They depend only
%     on the dip and are cached per dip (clear computeCornerflowVelocityfield
%     to reset). An error is raised if no solution meets the boundary
%     conditions.
%
% EXAMPLE
%   dip = pi/4;
%   [x, z] = meshgrid(linspace(0, 2, 50), linspace(-2, 0, 50));
%   [Vx, Vz] = computeCornerflowVelocityfield(dip, x(:), z(:), 1);
%
% SEE ALSO
%   getStrainratesLongterm, compute_cornerflow_velocityfield
%
% AUTHORS
%   Rishav Mallick (EOS), 2025

%% Input validation
if ~ismember(n, [1, 3])
    error('Power law exponent n must be 1 (Newtonian) or 3 (power-law)');
end

% Ensure column vectors
x = x(:);
z = z(:);

%% Compute polar coordinates
[theta, r2] = computePolarCoordinates(x, z);

%% Calculate velocity fields based on rheology
if n == 1
    % Newtonian fluid (analytical solution)
    [Vx, Vz] = computeNewtonianVelocity(dip, x, z, theta, r2);
else
    % Power-law fluid (n = 3, numerical solution)
    [Vx, Vz] = computePowerLawVelocity(dip, x, z, theta);
end

end

%% ========================================================================
%                         COORDINATE TRANSFORMATIONS
% =========================================================================

function [theta, r2] = computePolarCoordinates(x, z)
% Convert Cartesian to polar coordinates
% Note: atan2 returns angles in (-pi, pi]

theta = atan2(z, x);
r2 = x.^2 + z.^2;

% Avoid singularity at origin (r = 0)
r2(r2 == 0) = inf;

end

%% ========================================================================
%                      NEWTONIAN VELOCITY FIELD (n = 1)
% =========================================================================

function [Vx, Vz] = computeNewtonianVelocity(dip, x, z, theta, r2)
% Analytical solution for Newtonian corner flow

% Compute coefficients for both regions
[A_arc, B_arc, C_arc, D_arc] = getArcCoefficients(dip);
[A_ocn, B_ocn, C_ocn, D_ocn] = getOceanCoefficients(dip);

% Calculate velocities using stream function derivatives
vx_arc = calcNewtonianVx(A_arc, B_arc, C_arc, D_arc, x, z, theta, r2);
vz_arc = calcNewtonianVz(A_arc, B_arc, C_arc, D_arc, x, z, theta, r2);

vx_ocn = calcNewtonianVx(A_ocn, B_ocn, C_ocn, D_ocn, x, z, theta, r2);
vz_ocn = calcNewtonianVz(A_ocn, B_ocn, C_ocn, D_ocn, x, z, theta, r2);

% Combine results
Vx = [vx_arc, vx_ocn];
Vz = [vz_arc, vz_ocn];

end

function [A, B, C, D] = getArcCoefficients(dip)
% Continental arc (mantle wedge) coefficients
% Boundary conditions: Fixed top (v=0), Slab moves down-dip (v=U)

denom = dip^2 - sin(dip)^2;

C = (dip * sin(dip)) / denom;
D = (dip * cos(dip) - sin(dip)) / denom;

A = 0;        % Fixed top boundary
B = -C;       % Integration constant

end

function [A, B, C, D] = getOceanCoefficients(dip)
% Oceanic mantle (sub-slab) coefficients
% Boundary conditions: Plate moves in (v_r = -1), Slab moves out (v_r = 1)

alpha = pi - dip;  % Sub-slab corner angle

C = -sin(dip) / (alpha + sin(dip));
D = 2 * sin(dip/2)^2 / (alpha + sin(dip));

A = pi * C;
B = (dip - pi * cos(dip)) / (alpha + sin(dip));

end

function vx = calcNewtonianVx(A, B, C, D, x, z, theta, r2)
% Horizontal velocity component for Newtonian fluid

vx = -B - (C.*x.^2 + D.*x.*z) ./ r2 - D.*theta;

end

function vz = calcNewtonianVz(A, B, C, D, x, z, theta, r2)
% Vertical velocity component for Newtonian fluid

vz = A - (C.*x.*z + D.*z.^2) ./ r2 + C.*theta;

end

%% ========================================================================
%                    POWER-LAW VELOCITY FIELD (n = 3)
% =========================================================================

function [Vx, Vz] = computePowerLawVelocity(dip, x, z, theta)
% Power-law corner flow (n = 3). The coefficients depend only on the dip,
% so they are solved once per dip and cached (callers such as
% getStrainratesLongterm evaluate one triangle at a time).

persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
end
if isKey(cache, dip)
    coeff = cache(dip);
else
    % wedge: fixed upper plate on theta = 0, slab on theta = -dip
    coeff.arc = solvePowerLawCoefficients(0, [0, 0], -dip, [cos(dip), -sin(dip)], dip, 'mantle wedge');
    % oceanic mantle: incoming plate on theta = -pi, slab on theta = -dip
    coeff.ocn = solvePowerLawCoefficients(-pi, [1, 0], -dip, [cos(dip), -sin(dip)], dip, 'oceanic mantle');
    cache(dip) = coeff;
end

[vx_arc, vz_arc] = calcPowerLawVelocity(coeff.arc, theta, x, z);
[vx_ocn, vz_ocn] = calcPowerLawVelocity(coeff.ocn, theta, x, z);

Vx = [vx_arc, vx_ocn];
Vz = [vz_arc, vz_ocn];

end

function coeff = solvePowerLawCoefficients(t1, V1, t2, V2, dip, region)
% Coefficients [A, B, C, D] of the n = 3 solution for a corner between the
% rays theta = t1 and theta = t2 with prescribed velocities V1 and V2.
%
% The velocity is v(theta) = [-A, B] + C*w(theta; D) (independent of r),
% where w is the unit-amplitude solution with phase D. For a given D the
% boundary conditions are linear in A, B, C: subtracting them gives
%     C*(w(t2; D) - w(t1; D)) = V2 - V1,
% so D must make the two vectors parallel. That scalar equation is
% bracketed on a fine grid over one period of D and refined with fzero;
% C, A, B then follow directly. Near-singular roots (C large, the
% difference w(t2) - w(t1) small) are common for shallow dips, so every
% bracket is tried and only roots that meet the boundary conditions are
% kept. Distinct roots must give the same velocity field.

dV = V2 - V1;
period = 6*pi/sqrt(5);
Dgrid = linspace(0, period, 20001);
dw = unitSolution(t2, Dgrid) - unitSolution(t1, Dgrid);
crs = dw(:,1)*dV(2) - dw(:,2)*dV(1);
brackets = find(sign(crs(1:end-1)) ~= sign(crs(2:end)));

crossD = @(D) crossProduct(unitSolution(t2, D) - unitSolution(t1, D), dV);
tmid = (t1 + t2)/2;
sols = zeros(0, 4);
vmid = zeros(0, 2);
for k = brackets'
    D = fzero(crossD, Dgrid(k:k+1), optimset('TolX', eps));
    d = unitSolution(t2, D) - unitSolution(t1, D);
    C = (d*dV') / (d*d');
    if ~isfinite(C) || norm(C*d - dV) > 1e-8*norm(dV)
        continue    % sign change through a zero of w(t2) - w(t1), not a root
    end
    AB = V1 - C*unitSolution(t1, D);
    sols(end+1, :) = [-AB(1), AB(2), C, D]; %#ok<AGROW>
    vmid(end+1, :) = AB + C*unitSolution(tmid, D); %#ok<AGROW>
end

if isempty(sols)
    error('computeCornerflowVelocityfield:noSolution', ...
        'No n = 3 corner-flow solution meets the boundary conditions in the %s for dip = %.4g deg.', ...
        region, dip*180/pi);
end
if max(vecnorm(vmid - vmid(1,:), 2, 2)) > 1e-6
    error('computeCornerflowVelocityfield:ambiguous', ...
        'The n = 3 boundary conditions in the %s for dip = %.4g deg have more than one solution.', ...
        region, dip*180/pi);
end
coeff = sols(1, :);

end

function w = unitSolution(theta, D)
% w(theta; D): velocity [vx, vz] of the n = 3 solution with A = B = 0,
% C = 1 on the unit circle (it does not depend on r). Rows follow D.
D = D(:);
x = cos(theta);
z = sin(theta);
[vx, vz] = calcPowerLawVelocity([0, 0, 1, 0], theta + D, x, z);
w = [vx, vz];
end

function c = crossProduct(a, b)
c = a(1)*b(2) - a(2)*b(1);
end

function [vx, vz] = calcPowerLawVelocity(coeff, phi, x, z)
% Calculate velocity field for power-law fluid (n=3)
% Uses analytical form derived from stream function

A = coeff(1);
B = coeff(2);
C = coeff(3);
D = coeff(4);

% Compute radius
r = sqrt(x.^2 + z.^2);

% Angular terms
s5 = sqrt(5);
arg1 = s5/3 * (phi + D);
arg2 = s5 * (phi + D);

% Horizontal velocity component
vx = -A + C ./ r .* (...
    -27*z .* cos(arg1) + ...
      z .* cos(arg2) - ...
      s5*x .* (-9*sin(arg1) + sin(arg2)) ...
    );

% Vertical velocity component
vz = B - C ./ r .* (...
    -27*x .* cos(arg1) + ...
      x .* cos(arg2) + ...
      s5*z .* (-9*sin(arg1) + sin(arg2)) ...
    );

end