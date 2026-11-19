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
%   are found with fsolve from the boundary conditions.
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
%   - n = 3 needs the Optimization Toolbox (fsolve) and prints fsolve
%     diagnostics.
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
% Numerical solution for power-law corner flow (n = 3)

% Solve for coefficients using boundary conditions
coeff_arc = solveArcCoefficients(dip);
coeff_ocn = solveOceanCoefficients(dip);

% Calculate velocities
[vx_arc, vz_arc] = calcPowerLawVelocity(coeff_arc, theta, x, z);
[vx_ocn, vz_ocn] = calcPowerLawVelocity(coeff_ocn, theta, x, z);

% Combine results
Vx = [vx_arc, vx_ocn];
Vz = [vz_arc, vz_ocn];

end

function coeff = solveArcCoefficients(dip)
% Solve for mantle wedge coefficients [A, B, C, D] using fsolve
% Boundary conditions at phi=0 (surface) and phi=-dip (slab)

opts = optimoptions('fsolve', ...
    'Display', 'final-detailed', ...
    'MaxFunctionEvaluations', 1e4, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations',1e4);

x0 = [1, 1, 1, 0.1];  % Initial guess

coeff = fsolve(@residualArc, x0, opts);

    function R = residualArc(p)
        % Residual function for boundary conditions

        % Surface boundary (phi=0, z=0, x=1): velocity = 0
        [vx_surf, vz_surf] = calcPowerLawVelocity(p, 0, 1, 0);

        % Slab boundary (phi=-dip): velocity = [cos(dip), -sin(dip)]
        [vx_slab, vz_slab] = calcPowerLawVelocity(p, -dip, 1, -tan(dip));

        R = [vx_surf, ...
             vz_surf, ...
             vx_slab - cos(dip), ...
             vz_slab + sin(dip)];
    end

end

function coeff = solveOceanCoefficients(dip)
% Solve for oceanic mantle coefficients [A, B, C, D] using fsolve

opts = optimoptions('fsolve', ...
    'Display', 'final-detailed', ...
    'MaxFunctionEvaluations', 1e4, ...
    'FunctionTolerance', 1e-6,...
    'MaxIterations',1e4);

x0 = [1, 1, 1, 0.1];  % Initial guess

coeff = fsolve(@residualOcean, x0, opts);

    function R = residualOcean(p)
        % Residual function for boundary conditions

        % Incoming plate (phi=-pi, x=-1, z=0): velocity = [1, 0]
        [vx_plate, vz_plate] = calcPowerLawVelocity(p, -pi, -1, 0);

        % Slab boundary (phi=0, x=1, z=-tan(dip))
        [vx_slab, vz_slab] = calcPowerLawVelocity(p, -dip, 1, -tan(dip));

        R = [vx_plate - 1, ...
             vz_plate, ...
             vx_slab - cos(dip), ...
             vz_slab + sin(dip)];
    end

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