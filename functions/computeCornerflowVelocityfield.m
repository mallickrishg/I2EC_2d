function [Vx,Vz] = computeCornerflowVelocityfield(dip,x,z,n)
%COMPUTE_CORNERFLOW_VELOCITYFIELD Analytic corner flow for subduction zones.
%
%   [Vx,Vz] = COMPUTE_CORNERFLOW_VELOCITYFIELD(dip,x,z)
%
%   Computes the isoviscous corner flow velocity field for:
%   1. Continental Arc (Mantle Wedge): Fixed top, moving slab.
%   2. Oceanic Mantle (Sub-slab): Moving slab, incoming plate.
%
%   INPUTS:
%       dip - Slab dip angle in radians (e.g., pi/4)
%       x,z - Normalized coordinates [Nx1]. z is negative inside Earth.
%       n - Power law exponent (solutions only provided for: n = 1
%           Newtonian, n = 3)
%
%   OUTPUTS:
%       Vx, Vz - [Nx2] matrices. 
%                Column 1: Continental Arc (Wedge)
%                Column 2: Oceanic Mantle (Sub-slab)

% --- Setup Geometry ---
x = x(:);
z = z(:);

% Precompute polar coordinates for all points
% Note: atan2 returns values in (-pi, pi]
theta = atan2(z, x);
r2    = x.^2 + z.^2;

% Handle singularity at r=0 (trench) to avoid NaNs
r2(r2==0) = inf;

if n == 1
    % --- Continental Arc Coefficients ---
    % Domain: Wedge corner angle = dip
    % BCs: Fixed top (v=0), Slab moves down-dip (v=U)

    denom_arc = dip^2 - sin(dip)^2;

    % Analytical coefficients for Wedge (Standard Corner Flow)
    C_arc = (dip * sin(dip)) / denom_arc;
    D_arc = (dip * cos(dip) - sin(dip)) / denom_arc;

    % Integration constants for velocity form (A=0 for fixed top, B = -C)
    A_arc = 0;
    B_arc = -C_arc;

    % --- Oceanic Mantle Coefficients ---
    % Domain: Sub-slab corner angle beta = pi - dip
    % BCs: Plate moves in (v_r = -1), Slab moves out (v_r = 1)

    alpha = pi - dip;

    % Final Constants
    C_ocean = -sin(dip) / (alpha + sin(dip));
    A_ocean = pi * C_ocean;
    D_ocean = 2*(sin(dip/2))^2 / (alpha + sin(dip));
    B_ocean = (dip - pi*cos(dip)) / (alpha + sin(dip));

    % --- Compute Velocity Fields ---

    % Vectorized function for velocity components
    calc_vx = @(A,B,C,D) -B - (C.*x.^2 + D.*x.*z)./r2 - D.*theta;
    calc_vz = @(A,B,C,D)  A - (C.*x.*z + D.*z.^2)./r2 + C.*theta;

    % Evaluate Arc
    vx_arc = calc_vx(A_arc, B_arc, C_arc, D_arc);
    vz_arc = calc_vz(A_arc, B_arc, C_arc, D_arc);

    % Evaluate Ocean
    vx_ocn = calc_vx(A_ocean, B_ocean, C_ocean, D_ocean);
    vz_ocn = calc_vz(A_ocean, B_ocean, C_ocean, D_ocean);

    % --- Output ---
    Vx = [vx_arc, vx_ocn];
    Vz = [vz_arc, vz_ocn];

% For power-law rheology
elseif n == 3
    % --- Power-law (n = 3) Corner Flow ---        

    % === Continental arc coefficients ===
    % Unknowns: [A, B, C, D]

    coeff_arc = solve_coefficients_arc(dip);

    A_arc = coeff_arc(1);
    B_arc = coeff_arc(2);
    C_arc = coeff_arc(3);
    D_arc = coeff_arc(4);

    % === Oceanic mantle coefficients ===
    coeff_ocean = solve_coefficients_ocean(dip);

    A_ocean = coeff_ocean(1);
    B_ocean = coeff_ocean(2);
    C_ocean = coeff_ocean(3);
    D_ocean = coeff_ocean(4);
    

    % Evaluate Arc
    [vx_arc, vz_arc] = calc_powerlaw_velocity(A_arc,B_arc,C_arc,D_arc,theta,x,z);

    % Evaluate Ocean
    [vx_ocn, vz_ocn] = calc_powerlaw_velocity(A_ocean,B_ocean,C_ocean,D_ocean,theta,x,z);

    % --- Output ---
    Vx = [vx_arc, vx_ocn];
    Vz = [vz_arc, vz_ocn];

else
    error('Choose n = 1 or n = 3. Other analytical solutions are not yet implemented.')
end

end

%%%%%% Solve mantle wedge coefficients using fsolve %%%%%%%%
function coeff = solve_coefficients_arc(dip)

opts = optimoptions('fsolve','Display','final-detailed','MaxFunctionEvaluations',1e3,'FunctionTolerance',1e-6);

x0 = [1 1 1 0.1];

coeff = fsolve(@residual, x0, opts);

    function R = residual(p)
        A = p(1);
        B = p(2);
        C = p(3);
        D = p(4);

        % Arc surface (phi = 0, z = 0, x = 1)
        [vx1, vz1] = calc_powerlaw_velocity(A,B,C,D,0,1,0);

        % Slab surface (phi = -dip, z = -tan(dip), x = 1)
        [vx2, vz2] = calc_powerlaw_velocity(A,B,C,D,-dip,1,-tan(dip));

        R = [vx1, vz1, vx2 - cos(dip), vz2 + sin(dip)];
    end
end

%%%%%% Solve oceanic mantle coefficients using fsolve %%%%%%%%
function coeff = solve_coefficients_ocean(dip)

opts = optimoptions('fsolve','Display','final-detailed');

x0 = [1 1 1 0.1];

coeff = fsolve(@residual, x0, opts);

    function R = residual(p)
        A = p(1);
        B = p(2);
        C = p(3);
        D = p(4);

        % Oceanic slab surface (phi = -pi, z = 0, x = -1)
        [vx1, vz1] = calc_powerlaw_velocity(A,B,C,D,-pi,-1,0);

        % Slab surface (phi = -dip, z = -tan(dip), x = 1)
        [vx2, vz2] = calc_powerlaw_velocity(A,B,C,D,0,1,-tan(dip));
        % velocity_from_streamfunction(1, zslab, p);

        R = [vx1 - 1, vz1, vx2 - cos(dip), vz2 + sin(dip)];
    end
end

%%%%% Compute velocity-field for power-law fluid with known constants %%%%
function [vx,vz] = calc_powerlaw_velocity(A,B,C,D,phi,x,z)

vx = ...
    -A + ( ...
    -27*z.*cos(sqrt(5)/3*(phi + D)) ...
    +  z.*cos(sqrt(5)*(phi + D)) ...
    -  sqrt(5)*x.*( ...
    -9*sin(sqrt(5)/3*(phi + D)) ...
    +  sin(sqrt(5)*(phi + D)) ...
    ) ...
    ).*C./sqrt(x.^2 + z.^2);

vz = ...
    B - ( ...
    -27*x.*cos(sqrt(5)/3*(phi + D)) ...
    +  x.*cos(sqrt(5)*(phi + D)) ...
    +  sqrt(5)*z.*( ...
    -9*sin(sqrt(5)/3*(phi + D)) ...
    +  sin(sqrt(5)*(phi + D)) ...
    ) ...
    ).*C./sqrt(x.^2 + z.^2);


end