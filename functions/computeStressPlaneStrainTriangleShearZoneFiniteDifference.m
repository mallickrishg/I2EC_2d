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
%   - x3 must be non-negative (assertion).
%   - The closed-form displacement contains atan(N/D) terms whose
%     denominator vanishes on the lines through the edges of the triangle
%     and of its image above the surface; atan jumps by pi there although
%     the field is smooth (away from the triangle's own edges). A
%     finite-difference stencil across such a line gives a finite but
%     enormous strain (1e10 times too large), which the non-finite check
%     does not catch. Points within 1e-5 l of such a line (l = shortest
%     edge), and not on the triangle's own edges or near its vertices, are
%     therefore evaluated as the mean of two points shifted by +-1e-3 l
%     normal to the line (one-sided at the surface). This happens on
%     structured meshes aligned with the surface.
%
% SEE ALSO
%   computeStressPlaneStrainTriangleShearZone (closed form, preferred),
%   computeDisplacementPlaneStrainTriangleShearZone, computeShzStressKernelsBem
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

[s22,s23,s33]=stressDirect(x2,x3,A,B,C,e22,e23,e33,G,nu);

% points on the lines through the source or image edges (see NOTES)
X=[x2(:),x3(:)];
l=min([norm(A-B),norm(B-C),norm(C-A)]);
tol=1e-5*l;
delta=1e-3*l;
V=[A(:)';B(:)';C(:)'];
edges=[1 2;2 3;3 1];
flag=false(size(X,1),1);
nrm=zeros(size(X,1),2);
for img=[1,-1]
    P=V; P(:,2)=img*P(:,2);
    for k=1:3
        Q1=P(edges(k,1),:); Q2=P(edges(k,2),:);
        t=Q2-Q1; L=norm(t); n=[-t(2),t(1)]/L;
        d=(X-Q1)*n';
        u=(X-Q1)*t'/L^2;
        near=abs(d)<tol & ~flag ...
            & vecnorm(X-Q1,2,2)>2*delta & vecnorm(X-Q2,2,2)>2*delta;
        if img==1
            near=near & (u<0 | u>1);   % the triangle's own edge is a true discontinuity
        end
        flag(near)=true;
        nrm(near,:)=repmat(n,nnz(near),1);
    end
end
if any(flag)
    Xf=X(flag,:); nf=nrm(flag,:);
    Xp=Xf+delta*nf; Xm=Xf-delta*nf;
    % keep both points below the surface (one-sided where needed)
    up=Xp(:,2)<0; Xp(up,:)=Xm(up,:);
    dn=Xm(:,2)<0; Xm(dn,:)=Xp(dn,:);
    [a22,a23,a33]=stressDirect(Xp(:,1),Xp(:,2),A,B,C,e22,e23,e33,G,nu);
    [b22,b23,b33]=stressDirect(Xm(:,1),Xm(:,2),A,B,C,e22,e23,e33,G,nu);
    s22(flag)=(a22+b22)/2;
    s23(flag)=(a23+b23)/2;
    s33(flag)=(a33+b33)/2;
end

end

function [s22,s23,s33]=stressDirect(x2,x3,A,B,C,e22,e23,e33,G,nu)
% closed-form displacement differentiated by finite differences

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
