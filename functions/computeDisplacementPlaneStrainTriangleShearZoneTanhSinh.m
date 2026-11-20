function [u2,u3]=computeDisplacementPlaneStrainTriangleShearZoneTanhSinh( ...
    x2,x3,A,B,C,e22,e23,e33,nu,varargin)
% COMPUTEDISPLACEMENTPLANESTRAINTRIANGLESHEARZONETANHSINH  Displacement due to a strained triangle, by tanh-sinh quadrature.
%
%   [u2, u3] = computeDisplacementPlaneStrainTriangleShearZoneTanhSinh(x2, x3, A, B, C, e22, e23, e33, nu)
%   [u2, u3] = computeDisplacementPlaneStrainTriangleShearZoneTanhSinh(..., 'precision', h, 'bound', b)
%
%   Computes the same displacement as
%   computeDisplacementPlaneStrainTriangleShearZone by integrating the
%   half-space Green's functions along the three edges of the triangle
%   with double-exponential (tanh-sinh) quadrature. I2EC uses it only at
%   points where the closed form is not finite. Geometry (x3 is depth,
%   positive down):
%
%                 surface
%         -------------+-------------- E (x2)
%                      |
%                      |     + A
%                      |    /  .
%                      |   /     .
%                      |  /        .
%                      | /           .
%                      |/              + B
%                      /            .
%                     /|          /
%                    / :       .
%                   /  |    /
%                  /   : .
%                 /   /|
%                / .   :
%               +      |
%             C        :
%                      |
%                      D (x3)
%
% INPUTS
%   x2, x3          observation points: horizontal position and depth
%                   (m, x3 >= 0)
%   A, B, C         [1 x 2] triangle vertices [x2, x3] (m)
%   e22, e23, e33   eigenstrain components in the triangle (tensor
%                   components)
%   nu              Poisson's ratio
%
% OPTIONS
%   precision       quadrature step h  [0.001]
%   bound           range of the transformed variable; the sum runs over
%                   k = -n..n with n = fix(bound/h)  [3.5]
%
% OUTPUTS
%   u2              displacement in the x2 (east) direction, same size as
%                   x2 (m)
%   u3              displacement in the x3 (down) direction (m)
%
% NOTES
%   - The displacement depends on Poisson's ratio only, not on the shear
%     modulus.
%   - x3 must be non-negative (assertion).
%
% SEE ALSO
%   computeDisplacementPlaneStrainTriangleShearZone
%
% REFERENCES
%   Barbot, S. (2018). Deformation of a half-space from anelastic strain
%   confined in a tetrahedral volume. Bulletin of the Seismological Society
%   of America, 108(5A), 2687-2712. doi:10.1785/0120180058
%
% AUTHORS
%   Sylvain Barbot (sbarbot@ntu.edu.sg) - Feb 9, 2018, Los Angeles.

assert(min(x3(:))>=0,'depth must be positive.');

% process optional input
p = inputParser;
p.addParameter('precision',0.001,@validatePrecision);
p.addParameter('bound',3.5,@validateBound);
p.parse(varargin{:});
optionStruct = p.Results;

% Lame parameter
lambda=2*nu/(1-2*nu);

% array size
s=size(x2);

% isotropic strain
ekk=e22+e33;

% unit vectors
nA = [C(2)-B(2);
      B(1)-C(1)]/norm(C-B);
nB = [C(2)-A(2);
      A(1)-C(1)]/norm(C-A);
nC = [B(2)-A(2);
      A(1)-B(1)]/norm(B-A);
  
% check that unit vectors are pointing outward
if (nA'*(A(:)-(B(:)+C(:))/2))>0
    nA=-nA;
end
if (nB'*(B(:)-(A(:)+C(:))/2))>0
    nB=-nB;
end
if (nC'*(C(:)-(A(:)+B(:))/2))>0
    nC=-nC;
end

% Radii
r1=@(y2,y3) sqrt((x2-y2).^2+(x3-y3).^2);
r2=@(y2,y3) sqrt((x2-y2).^2+(x3+y3).^2);

% parameterized line integral
y2=@(t,A,B) (A(1)+B(1))/2+t*(B(1)-A(1))/2;
y3=@(t,A,B) (A(2)+B(2))/2+t*(B(2)-A(2))/2;

% Green's functions
G22=@(y2,y3) -1/(2*pi*(1-nu))*( ...
    (3-4*nu)/4*log(r1(y2,y3))+(8*nu^2-12*nu+5)/4*log(r2(y2,y3))+(x3-y3).^2./r1(y2,y3).^2/4 ...
    +((3-4*nu)*(x3+y3).^2+2*y3.*(x3+y3)-2*y3.^2)./r2(y2,y3).^2/4-(y3.*x3.*(x3+y3).^2)./r2(y2,y3).^4 ...
    );
G23=@(y2,y3) 1/(2*pi*(1-nu))*( ...
    (1-2*nu)*(1-nu)*atan((x2-y2)./(x3+y3))+((x3-y3).*(x2-y2))./r1(y2,y3).^2/4 ...
    +(3-4*nu)*((x2-y2).*(x3-y3))./r2(y2,y3).^2/4-(y3.*x3.*(x2-y2).*(x3+y3))./r2(y2,y3).^4 ...
    );

G32=@(y2,y3) 1/(2*pi*(1-nu))*( ...
    -(1-2*nu)*(1-nu)*atan((x2-y2)./(x3+y3)) ...
    +(x3-y3).*(x2-y2)./r1(y2,y3).^2/4 ...
    +(3-4*nu)*(x2-y2).*(x3-y3)./r2(y2,y3).^2/4 ...
    +y3.*x3.*(x2-y2).*(x3+y3)./r2(y2,y3).^4 ...
    );
G33=@(y2,y3) 1/(2*pi*(1-nu))*( ...
    -(3-4*nu)/4*log(r1(y2,y3))-(8*nu^2-12*nu+5)/4*log(r2(y2,y3)) ...
    -(x2-y2).^2./r1(y2,y3).^2/4+(2*y3.*x3-(3-4*nu)*(x2-y2).^2)./r2(y2,y3).^2/4-y3.*x3.*(x2-y2).^2./r2(y2,y3).^4 ...
     );

% moment density
m22=lambda*ekk+2*e22;
m23=2*e23;
m33=lambda*ekk+2*e33;
 
% function IU2 is the integrand for displacement component u2
IU2=@(t) ...
    (m22*nC(1)+m23*nC(2))*norm(B-A)/2*G22(y2(t,A,B),y3(t,A,B)) ...
   +(m23*nC(1)+m33*nC(2))*norm(B-A)/2*G32(y2(t,A,B),y3(t,A,B)) ...
   +(m22*nA(1)+m23*nA(2))*norm(C-B)/2*G22(y2(t,B,C),y3(t,B,C)) ...
   +(m23*nA(1)+m33*nA(2))*norm(C-B)/2*G32(y2(t,B,C),y3(t,B,C)) ...
   +(m22*nB(1)+m23*nB(2))*norm(A-C)/2*G22(y2(t,C,A),y3(t,C,A)) ...
   +(m23*nB(1)+m33*nB(2))*norm(A-C)/2*G32(y2(t,C,A),y3(t,C,A));

% function IU3 is the integrand for displacement component u3
IU3=@(t) ...
    (m22*nC(1)+m23*nC(2))*norm(B-A)/2*G23(y2(t,A,B),y3(t,A,B)) ...
   +(m23*nC(1)+m33*nC(2))*norm(B-A)/2*G33(y2(t,A,B),y3(t,A,B)) ...
   +(m22*nA(1)+m23*nA(2))*norm(C-B)/2*G23(y2(t,B,C),y3(t,B,C)) ...
   +(m23*nA(1)+m33*nA(2))*norm(C-B)/2*G33(y2(t,B,C),y3(t,B,C)) ...
   +(m22*nB(1)+m23*nB(2))*norm(A-C)/2*G23(y2(t,C,A),y3(t,C,A)) ...
   +(m23*nB(1)+m33*nB(2))*norm(A-C)/2*G33(y2(t,C,A),y3(t,C,A));

% numerical solution
h=optionStruct.precision;
n=fix(1/h*optionStruct.bound);
u2=zeros(s);
u3=zeros(s);
for k=-n:n
    wk=(0.5*h*pi*cosh(k*h))./(cosh(0.5*pi*sinh(k*h))).^2;
    xk=tanh(0.5*pi*sinh(k*h));
    u2=u2+wk*IU2(xk);
    u3=u3+wk*IU3(xk);
end

end

function p = validatePrecision(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error(message('MATLAB:mcmc:invalidPrecision'));
end
p = true;
end

function p = validateBound(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error(message('MATLAB:mcmc:invalidBound'));
end
p = true;
end
