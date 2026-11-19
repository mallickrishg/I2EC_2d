function [Stress] = LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E)
% LDSTRESSHS  Half-space stress from a line dislocation (crack).
%
%   Stress = geometry.LDstressHS(x, y, xe, ye, a, Beta, Ds, Dn, nu, E)
%
%   Stress at the observation points due to a single straight element of
%   a 2-D displacement-discontinuity boundary element, with uniform shear
%   and opening, in a homogeneous elastic half-space whose free surface is
%   y = 0 (source plus image terms).
%
% INPUTS
%   x, y    [Nobs x 1] observation points (m); y must be <= 0
%   xe, ye  [1 x 1] element midpoint (m); ye must be below 0
%   a       [1 x 1] element half-length (m)
%   Beta    [1 x 1] angle of the element from the x axis (radians); 0 when
%           the normal points down (-y). I2EC passes -deg2rad(dip).
%   Ds      [1 x 1] shear displacement discontinuity (m); positive is left
%           lateral
%   Dn      [1 x 1] normal displacement discontinuity (m); positive is
%           opening
%   nu      [1 x 1] Poisson's ratio
%   E       [1 x 1] Young's modulus (MPa); I2EC passes 2*G*(1 + nu)
%
% OUTPUTS
%   Stress  [Nobs x 3] stress [Sxx, Syy, Sxy] at the observation points
%           (same units as E)
%
% NOTES
%   - Stops with an error if any observation point has y > 0.
%
% EXAMPLE
%   [x, y] = meshgrid(-2:0.05:2, -4:0.05:0);
%   Stress = geometry.LDstressHS(x(:), y(:), 0, -2, 1, 0, 0, 1, 0.25, 1);
%   scatter(x(:), y(:), 10, Stress(:,3), 'filled'); colorbar
%
% SEE ALSO
%   geometry.LDdispHS, geometry.computeFaultTractionKernels,
%   geometry.computeFaultStressKernels
%
% REFERENCES
%   Crouch, S. L. & Starfield, A. M. (1983). Boundary Element Methods in
%   Solid Mechanics.
%   Ritz, E. & Pollard, D. D. (2012) (notation).
%
% AUTHORS
%   Tim Davis/Steve Martel
%   Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%   Modified form of Steve Martels BEM scripts from:
%   http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%   Also includes some notation from Ritz & Pollard 2012.
%   Added to I2EC by Rishav Mallick, 2023


if any(y>0)
    error('Half-space solution: Z coordinates must be negative!')
end
% The shear modulus, sm, is related to the prescribed elastic constants.
sm = E/(2*(1+nu));
% Define material constant used in calculating influence coefficients.
con = 1/(4*pi*(1-nu)); 
cons = 2*sm;
Dxb = Ds; Dyb = -Dn;
sb = sin(Beta); cb = cos(Beta);
s2b = sin(2*Beta); c2b = cos(2*Beta);
s3b = sin(3*Beta); c3b = cos(3*Beta);
s4b = sin(4*Beta); c4b = cos(4*Beta);

% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
XB = (x-xe)*cb + (y-ye)*sb;
YB = -(x-xe)*sb + (y-ye)*cb;

% Coordinates of the image dislocation
XBi = (x-xe)*cb - (y+ye)*sb;		%equation 7.4.6 C&S
YBi = (x-xe)*sb + (y+ye)*cb ;

% Fix roundoff errors in Ybi and Yb from trig function problems
bad = abs(YBi)<1e-10; %Steve Martels Fix
bad2 = abs(YB)<1e-10;
YBi(bad) = 0;
YB(bad2) = 0;
clear bad bad2


% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%   which are used to calculate the displacement and stress components. 
% It is understood that X and Y refer to XB and YB.
% First abbreviate repeated terms in the derivatives of f(x,y):
Y2 = YB.^2;
XMa = XB-a; XPa = XB+a; 
XMa2 = XMa.^2; XPa2 = XPa.^2; 
R1S = XMa2 + Y2; R1S2 = R1S.^2; 
R2S = XPa2 + Y2; R2S2 = R2S.^2;

% Same thing for the image dislocation
Y2i = YBi.^2;
XMai = XBi-a; XPai = XBi+a; 
XMa2i = XMai.^2; XPa2i = XPai.^2; 
R1Si = XMa2i + Y2i; R1S2i = R1Si.^2; 
R2Si = XPa2i + Y2i; R2S2i = R2Si.^2;

FF4 = con*(YB./R1S - YB./R2S); 
FF5 = con*(XMa./R1S - XPa./R2S);
% The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
FF6 = con*((XMa2 - Y2)./R1S2 - (XPa2 - Y2)./R2S2);
FF7 = 2*con*YB.*(XMa./R1S2 - XPa./R2S2);


FF4i = con*(YBi./R1Si - YBi./R2Si); 
FF5i = con*(XMai./R1Si - XPai./R2Si);


% The halfspace examples of eqs. 5.5.3a and b of C&S, p. 91.
% See Appendix A of: Martel, S.J. and Langley, J.S., 2006. Propagation of
% normal faults to the surface in basalt, Koae fault system, Hawaii.
% Journal of Structural Geology, 28(12), pp.2123-2143.
FF6i = con*((XMa2i - Y2i)./R1S2i - (XPa2i - Y2i)./R2S2i);
FF7i = 2*con*YBi.*(XMai./R1S2i - XPai./R2S2i);

%*Tim* I used MATLABs symbolic to find these not eq's A.3 and A.4 of Martel
% Used Eq.A.1 on variable FF7i (expanded).
FF8i =(YBi.*(1./((a + XBi).^2 + YBi.^2).^2 - 1./(YBi.^2 + (a - XBi).^2).^2 + (2.*(a - XBi).*(2.*a - 2.*XBi))./(YBi.^2 + (a - XBi).^2).^3 - (2.*(a + XBi).*(2.*a + 2.*XBi))./((a + XBi).^2 + YBi.^2).^3))./(2.*pi.*(nu - 1));
FF9i =((a - XBi)./(YBi.^2 + (a - XBi).^2).^2 + (a + XBi)./((a + XBi).^2 + YBi.^2).^2)./(2.*pi.*(nu - 1)) - (YBi.*((4.*YBi.*(a + XBi))./((a + XBi).^2 + YBi.^2).^3 + (4.*YBi.*(a - XBi))./(YBi.^2 + (a - XBi).^2).^3))./(2.*pi.*(nu - 1));



% Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
Sxx =  cons*Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
 +cons*Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
Syy =  cons*Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
 +cons*Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
Sxy =  cons*Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
 +cons*Dyb*(-YB.*(c2b*FF6 - s2b*FF7)); 

%  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit SHEAR and
%  NORMAL displacement discontinuity
Sxxi_s = cons*Dxb*(FF4i - 3.*(c2b.*FF4i - s2b.*FF5i) +....
(2.*y.*(cb - 3.*c3b) + 3.*YB.*c2b).*FF6i +.....
(2.*y.*(sb - 3.*s3b) + 3.*YB.*s2b).*FF7i -.....
2.*y.*(y.*c4b - YB.*c3b).*FF8i -....
2.*y.*(y.*s4b - YB.*s3b).*FF9i);

Sxxi_n = cons*Dyb*(FF5i + (2.*y.*(sb - 2.*s3b) +.... 
3.*YB.*s2b).*FF6i - (2.*y.*(cb - 2.*c3b) +....
3.*YB.*c2b).*FF7i - 2.*y.*(y.*s4b - YB.*s3b).*FF8i +....
2.*y.*(y.*c4b - YB.*c3b).*FF9i);

Syyi_s = cons*Dxb*(FF4i - (c2b.*FF4i - s2b.*FF5i) -....
(4.*y.*sb.*s2b - YB.*c2b).*FF6i +.....
(4.*y.*sb.*c2b + YB.*s2b).*FF7i +....
2.*y.*(y.*c4b - YB.*c3b).*FF8i +....
2.*y.*(y.*s4b - YB.*s3b).*FF9i);

Syyi_n = cons*Dyb*(FF5i - (2.*y.*sb - YB.*s2b).*FF6i +.... 
(2.*y.*cb - YB.*c2b).*FF7i +....
2.*y.*(y.*s4b - YB.*s3b).*FF8i -.... 
2.*y.*(y.*c4b - YB.*c3b).*FF9i);

Sxyi_s = cons*Dxb*(s2b.*FF4i + c2b.*FF5i +.... 
(2.*y.*sb.*(1+4.*c2b) - YB.*s2b).*FF6i +.....
(2.*y.*cb.*(3-4.*c2b) + YB.*c2b).*FF7i +..... 
2.*y.*(y.*s4b - YB.*s3b).*FF8i -.....
2.*y.*(y.*c4b - YB.*c3b).*FF9i);

Sxyi_n = cons*Dyb*((4.*y.*sb.*s2b + YB.*c2b).*FF6i -.... 
(4.*y.*sb.*c2b - YB.*s2b).*FF7i -.....
2.*y.*(y.*c4b - YB.*c3b).*FF8i -.... 
2.*y.*(y.*s4b - YB.*s3b).*FF9i); 


Sxxi=Sxxi_s+Sxxi_n;
Syyi=Syyi_s+Syyi_n;
Sxyi=Sxyi_s+Sxyi_n;

Sxx=Sxx+Sxxi;
Syy=Syy+Syyi;
Sxy=Sxy+Sxyi;


Stress=[Sxx(:),Syy(:),Sxy(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
