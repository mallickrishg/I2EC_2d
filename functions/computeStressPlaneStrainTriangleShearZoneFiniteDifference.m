function [s22,s23,s33]=computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
    x2,x3,A,B,C,e22,e23,e33,G,nu)
% function COMPUTESTRESSPLANESTRAINTRIANGLESHEARZONE computes the
% stress field associated with deforming triangle strain volume
% considering the following geometry using the analytic solution of
% Barbot (2018).
%
%              surface
%      -------------+-------------- E (x2)
%                   |
%                   |     + A
%                   |    /  . 
%                   |   /     .  
%                   |  /        .            
%                   | /           .      
%                   |/              + B
%                   /            .
%                  /|          /  
%                 / :       .
%                /  |    /
%               /   : .
%              /   /|
%             / .   :
%            +      |
%          C        :
%                   |
%                   D (x3)
%
%
% Input:
% x2, x3             east coordinates and depth of the observation point,
% A, B, C            east and depth coordinates of the vertices,
% eij                source strain component 22, 23 and 33 in the strain
%                    volume,
% G, nu              rigidity and Poisson's ratio in the half space.
%
% Output:
% s22                horizontal uniaxial stress,
% s23                shear stress
% s33                vertical uniaxial stress.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - March 7, 2018, Singapore.

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
