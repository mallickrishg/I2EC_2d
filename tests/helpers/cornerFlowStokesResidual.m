function [stokes, radial] = cornerFlowStokesResidual(velfun, theta, n)
%CORNERFLOWSTOKESRESIDUAL  Check that a corner flow solves power-law Stokes.
%
%   [stokes, radial] = cornerFlowStokesResidual(velfun, theta, n)
%
%   velfun(x, z) returns [vx, vz] (column vectors); theta is a vector of
%   angles spanning one region of the corner, away from its boundaries.
%
%   A corner flow has a velocity that depends on the angle only, v(theta)
%   (stream function psi = r g(theta), v_r = g', v_theta = -g). Its only
%   strain rate is edot_rtheta = F/(2r) with F = v_r' - v_theta (= g'' + g).
%   With deviatoric stress tau = 2 eta edot and eta = edot_II^((1-n)/n),
%   tau_rtheta = r^(-1/n) h(theta) with h = sign(F)|F|^(1/n), and the two
%   momentum equations reduce (after eliminating the pressure) to
%       h'' + k^2 h = 0,   k^2 = (2n - 1)/n^2,
%   i.e. k = 1 for n = 1 and k = sqrt(5)/3 for n = 3. So the flow solves
%   Stokes iff h is a sinusoid of wavenumber k.
%
%   stokes  max deviation of h from its least-squares fit by
%           cos(k theta), sin(k theta), relative to max|h|, over the points
%           where |F| > 1e-3 max|F| (the cube root amplifies rounding
%           where F = 0). ~1e-5 or less for a solution, O(0.1-1) otherwise.
%   radial  max change of the velocity between r = 1 and r = 7.3 (it must
%           not depend on r).
%
%   F uses one central difference in theta (step 1e-5), on the unit circle.

theta = theta(:);
dth = 1e-5;
k = sqrt(2*n - 1) / n;

[vr, vt] = polarVelocity(velfun, theta, 1);
vrp = polarVelocity(velfun, theta + dth, 1);
vrm = polarVelocity(velfun, theta - dth, 1);
F = (vrp - vrm) / (2*dth) - vt;
h = sign(F) .* abs(F).^(1/n);

M = [cos(k*theta), sin(k*theta)];
c = M \ h;
keep = abs(F) > 1e-3 * max(abs(F));
stokes = max(abs(M(keep,:)*c - h(keep))) / max(abs(h));

[vr2, vt2] = polarVelocity(velfun, theta, 7.3);
radial = max(abs([vr2 - vr; vt2 - vt]));
end

function [vr, vt] = polarVelocity(velfun, theta, r)
[vx, vz] = velfun(r*cos(theta), r*sin(theta));
vr = vx.*cos(theta) + vz.*sin(theta);
vt = -vx.*sin(theta) + vz.*cos(theta);
end
