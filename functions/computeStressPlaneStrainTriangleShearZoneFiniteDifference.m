function [s22,s23,s33]=computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
    x2,x3,A,B,C,e22,e23,e33,G,nu)
% COMPUTESTRESSPLANESTRAINTRIANGLESHEARZONEFINITEDIFFERENCE  Stress due to a strained triangle, by finite differences.
%
%   [s22, s23, s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference(x2, x3, A, B, C, e22, e23, e33, G, nu)
%
%   Computes the plane-strain stress at (x2, x3) due to uniform
%   eigenstrain (e22, e23, e33) in the triangle ABC in a homogeneous
%   elastic half-space. The strain is obtained by numerically
%   differentiating the analytic displacement of Barbot (2018)
%   (computeDisplacementPlaneStrainTriangleShearZone); the eigenstrain is
%   removed inside the triangle and Hooke's law gives the stress. This is
%   the stress used for the I2EC shear-zone kernels (through
%   geometry.computeShzStressKernels). Geometry (x3 is depth, positive
%   down):
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
%   G               shear modulus (MPa)
%   nu              Poisson's ratio
%
% OUTPUTS
%   s22             horizontal normal stress, same size as x2 (MPa)
%   s23             shear stress (MPa)
%   s33             vertical normal stress (MPa)
%
% NOTES
%   - Step h = 1e-6 times the shortest edge of the triangle: central
%     differences in x2, a forward difference (u(x3 + h) - u(x3))/h in x3.
%   - The eigenstrain is subtracted at points strictly inside the
%     triangle; I2EC evaluates the stress at cell centroids.
%   - Where the closed-form displacement is not finite, it is recomputed
%     by tanh-sinh quadrature (see computeDisplacementPlaneStrainTriangleShearZone).
%   - Used instead of the closed-form stress of
%     computeStressPlaneStrainTriangleShearZone.
%   - x3 must be non-negative (assertion).
%
% SEE ALSO
%   computeDisplacementPlaneStrainTriangleShearZone,
%   computeStressPlaneStrainTriangleShearZone, computeShzStressKernelsBem
%
% REFERENCES
%   Barbot, S. (2018). Deformation of a half-space from anelastic strain
%   confined in a tetrahedral volume. Bulletin of the Seismological Society
%   of America, 108(5A), 2687-2712. doi:10.1785/0120180058
%
% AUTHORS
%   Sylvain Barbot (sbarbot@ntu.edu.sg) - March 7, 2018, Singapore
%   Added to I2EC by Rishav Mallick, 2023

assert(min(x3(:))>=0,'depth must be positive.');

% Lame parameter
lambda=2*G*nu/(1-2*nu);

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

l=min([sqrt((A(1)-B(1))^2+(A(2)-B(2))^2), ...
       sqrt((A(1)-C(1))^2+(A(2)-C(2))^2), ...
       sqrt((B(1)-C(1))^2+(B(2)-C(2))^2), ...
    ]);

step=1e-6;

p2=l*step;
p3=l*step;

[u2p2,u3p2]=computeDisplacementPlaneStrainTriangleShearZone( ...
    x2+p2,x3,A,B,C,e22,e23,e33,nu);
[u2m2,u3m2]=computeDisplacementPlaneStrainTriangleShearZone( ...
    x2-p2,x3,A,B,C,e22,e23,e33,nu);

[u2p3,u3p3]=computeDisplacementPlaneStrainTriangleShearZone( ...
    x2,x3+p3,A,B,C,e22,e23,e33,nu);
[u2m3,u3m3]=computeDisplacementPlaneStrainTriangleShearZone( ...
    x2,x3   ,A,B,C,e22,e23,e33,nu);

% strain components
E22=(u2p2-u2m2)/(2*l*step);
E23=((u2p3-u2m3)/(l*step)+(u3p2-u3m2)/(2*l*step))/2;
E33=(u3p3-u3m3)/(l*step);

% remove anelastic strain
Omega=@(x2,x3) heaviside(((A(1)+B(1))/2-x2)*nC(1)+((A(2)+B(2))/2-x3)*nC(2)) ...
             .*heaviside(((B(1)+C(1))/2-x2)*nA(1)+((B(2)+C(2))/2-x3)*nA(2)) ...
             .*heaviside(((C(1)+A(1))/2-x2)*nB(1)+((C(2)+A(2))/2-x3)*nB(2));

E22=E22-Omega(x2,x3)*e22;
E23=E23-Omega(x2,x3)*e23;
E33=E33-Omega(x2,x3)*e33;

% stress components
s22=lambda*(E22+E33)+2*G*E22;
s23=2*G*E23;
s33=lambda*(E22+E33)+2*G*E33;

    function y=heaviside(x)
        y=x>0;
    end
end
