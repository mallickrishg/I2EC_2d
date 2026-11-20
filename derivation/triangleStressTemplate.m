function [s22,s23,s33]=computeStressPlaneStrainTriangleShearZone( ...
    x2,x3,A,B,C,e22,e23,e33,G,nu)
% COMPUTESTRESSPLANESTRAINTRIANGLESHEARZONE  Stress due to a strained triangle, in closed form.
%
%   [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(x2, x3, A, B, C, e22, e23, e33, G, nu)
%
%   Computes the plane-strain stress at (x2, x3) due to uniform
%   eigenstrain (e22, e23, e33) in the triangle ABC in a homogeneous
%   elastic half-space. The displacement gradient is the exact derivative
%   of the closed-form displacement of Barbot (2018)
%   (computeDisplacementPlaneStrainTriangleShearZone), generated with the
%   Symbolic Math Toolbox by derivation/generateTriangleStress.m; the
%   eigenstrain is removed inside the triangle and Hooke's law gives the
%   stress. It replaces an earlier finite-difference version (removed
%   2026-11-20; truncation errors ~1e-6, up to 1e-2 near edge lines).
%   Geometry (x3 is depth, positive down):
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
%                   components): scalars, or [1 x K] rows for K source
%                   strains at once (the edge integrals are shared)
%   G               shear modulus (MPa)
%   nu              Poisson's ratio
%
% OUTPUTS
%   s22             horizontal normal stress, same size as x2 (MPa); for
%                   K source strains, [numel(x2) x K]
%   s23             shear stress (MPa)
%   s33             vertical normal stress (MPa)
%
% NOTES
%   - d/dx atan(N/D) is written (N'D - ND')/(N^2 + D^2), so the derivative
%     is smooth across D = 0, where the atan itself jumps by pi.
%   - Removable singularities. Some terms of the closed form behave like
%     atan(N/D)/D and cancel only between the two ends of an edge, on the
%     line through the edge of the image triangle and, for horizontal or
%     vertical edges, on the line through the edge itself. Within a
%     distance d of such a line the closed form loses accuracy like
%     1e-13*(r/d)^2 (r = distance from the edge, at least its length L)
%     and is NaN on it. There, that edge's contribution is instead
%     integrated numerically: quadrature of the differentiated Green's
%     function along the edge (64-point Gauss-Legendre farther than L/4
%     from the edge, tanh-sinh nearer), for points with d < 1e-2 r.
%   - Points where the result is still not finite (on the triangle's
%     edges and vertices, where the stress is singular or discontinuous)
%     get the quadrature value of every edge.
%   - The eigenstrain is subtracted at points strictly inside the triangle.
%   - Edges within 1e-7 of horizontal or vertical are snapped, as in
%     computeDisplacementPlaneStrainTriangleShearZone.
%   - x3 must be non-negative (assertion).
%
% SEE ALSO
%   computeDisplacementPlaneStrainTriangleShearZone,
%   computeDisplacementPlaneStrainTriangleShearZoneTanhSinh,
%   geometry.computeShzStressKernels
%
% REFERENCES
%   Barbot, S. (2018). Deformation of a half-space from anelastic strain
%   confined in a tetrahedral volume. Bulletin of the Seismological Society
%   of America, 108(5A), 2687-2712. doi:10.1785/0120180058
%
% AUTHORS
%   Closed-form displacement: Sylvain Barbot (sbarbot@ntu.edu.sg), 2018.
%   Analytic derivative: Rishav Mallick, 2026

assert(min(x3(:))>=0,'depth must be positive.');

sz=size(x2);
x2=x2(:); x3=x3(:);

% make nearly horizontal / vertical edges exactly so
V=snapNearlyAlignedEdges([A(:)';B(:)';C(:)']);
A=V(1,:); B=V(2,:); C=V(3,:);

% Lame parameter (per unit shear modulus) and moment density, one
% column per source strain
e22=e22(:)'; e23=e23(:)'; e33=e33(:)';
lambda=2*nu/(1-2*nu);
ekk=e22+e33;
m22=lambda*ekk+2*e22;
m23=2*e23;
m33=lambda*ekk+2*e33;

% outward unit normals
nA=[C(2)-B(2); B(1)-C(1)]/norm(C-B);
nB=[C(2)-A(2); A(1)-C(1)]/norm(C-A);
nC=[B(2)-A(2); A(1)-B(1)]/norm(B-A);
if (nA'*(A(:)-(B(:)+C(:))/2))>0, nA=-nA; end
if (nB'*(B(:)-(A(:)+C(:))/2))>0, nB=-nB; end
if (nC'*(C(:)-(A(:)+B(:))/2))>0, nC=-nC; end

% displacement gradient due to each edge, per unit edge traction; the
% columns are the x2- and x3-derivatives of T22, T32, T23, T33, with
% u2 = c2*T22 + c3*T32 and u3 = c2*T23 + c3*T33
edges={A,B,nC; B,C,nA; C,A,nB};
base=cell(3,1);
for k=1:3
    P=edges{k,1}; Q=edges{k,2};
    base{k}=edgeGradient(x2,x3,P,Q,nu);
    % points near a line where the closed form cancels (see NOTES)
    bad=nearRemovableSingularity(x2,x3,P,Q) | any(~isfinite(base{k}),2);
    if any(bad)
        base{k}(bad,:)=edgeGradientQuadrature(x2(bad),x3(bad),P,Q,nu);
    end
end
% on the triangle's own edges and vertices: quadrature for every edge
bad=any(~isfinite(base{1}+base{2}+base{3}),2);
if any(bad)
    for k=1:3
        base{k}(bad,:)=edgeGradientQuadrature(x2(bad),x3(bad),edges{k,1},edges{k,2},nu);
    end
end

% displacement gradient for each source strain [N x K]
u22=0; u23=0; u32=0; u33=0;
for k=1:3
    n=edges{k,3}; b=base{k};
    c2=m22*n(1)+m23*n(2);       % traction on the edge, x2 component [1 x K]
    c3=m23*n(1)+m33*n(2);       % x3 component
    u22=u22+b(:,1)*c2+b(:,3)*c3;
    u23=u23+b(:,2)*c2+b(:,4)*c3;
    u32=u32+b(:,5)*c2+b(:,7)*c3;
    u33=u33+b(:,6)*c2+b(:,8)*c3;
end

% strain, minus the eigenstrain inside the triangle
in=insideTriangle(x2,x3,A,B,C);
E22=u22-in*e22;
E23=(u23+u32)/2-in*e23;
E33=u33-in*e33;

% Hooke's law
if 1==numel(e22), osz=sz; else, osz=[numel(x2),numel(e22)]; end
s22=reshape(G*(lambda*(E22+E33)+2*E22),osz);
s23=reshape(G*2*E23,osz);
s33=reshape(G*(lambda*(E22+E33)+2*E33),osz);

end

function b=edgeGradient(x2,x3,P,Q,nu)
% closed-form x2- and x3-derivatives of T22, T32, T23, T33 for edge PQ,
% T = I(h) - I(-h)
a=(Q-P)/norm(Q-P); m=(P+Q)/2; h=norm(Q-P)/2;
a2=a(1); a3=a(2); m2=m(1); m3=m(2);
if 0==a2
    f22=@dI22_v; f33=@dI33_v; f23=@dI23_v; f32=@dI32_v;
elseif 0==a3
    f22=@dI22_g; f33=@dI33_g; f23=@dI23_h; f32=@dI32_h;
else
    f22=@dI22_g; f33=@dI33_g; f23=@dI23_g; f32=@dI32_g;
end
[p22,q22]=f22(x2,x3,m2,m3,a2,a3,h,nu); [r22,s22]=f22(x2,x3,m2,m3,a2,a3,-h,nu);
[p32,q32]=f32(x2,x3,m2,m3,a2,a3,h,nu); [r32,s32]=f32(x2,x3,m2,m3,a2,a3,-h,nu);
[p23,q23]=f23(x2,x3,m2,m3,a2,a3,h,nu); [r23,s23]=f23(x2,x3,m2,m3,a2,a3,-h,nu);
[p33,q33]=f33(x2,x3,m2,m3,a2,a3,h,nu); [r33,s33]=f33(x2,x3,m2,m3,a2,a3,-h,nu);
b=[p22-r22, q22-s22, p32-r32, q32-s32, p23-r23, q23-s23, p33-r33, q33-s33];
end

function bad=nearRemovableSingularity(x2,x3,P,Q)
% points within the angle 1e-2 (seen from the edge midpoint, and at
% least 1e-2 L wide) of the image edge line and, for axis-aligned edges,
% of the edge line itself (L = edge length)
tol=1e-2;
L=norm(Q-P); a=(Q-P)/L; n=[-a(2),a(1)];
w=tol*max(L,hypot(x2-(P(1)+Q(1))/2,x3-(P(2)+Q(2))/2));
% image edge (x3 -> -x3)
nI=[-a(2),-a(1)]; PI=[P(1),-P(2)];
bad=abs((x2-PI(1))*nI(1)+(x3-PI(2))*nI(2))<w;
if 0==a(1) || 0==a(2)
    bad=bad | abs((x2-P(1))*n(1)+(x3-P(2))*n(2))<w;
end
end

function d=edgeGradientQuadrature(x2,x3,P,Q,nu)
% the columns of edgeGradient by quadrature of the differentiated Green's
% function: 64-point Gauss-Legendre for points farther than L/4 from the
% edge, tanh-sinh (step 5e-4) for nearer ones
L=norm(Q-P);
tt=max(0,min(1,((x2-P(1))*(Q(1)-P(1))+(x3-P(2))*(Q(2)-P(2)))/L^2));
far=hypot(x2-P(1)-tt*(Q(1)-P(1)),x3-P(2)-tt*(Q(2)-P(2)))>L/4;
d=zeros(numel(x2),8);
[t,w]=gaussLegendre64();
d(far,:)=edgeQuadrature(x2(far),x3(far),P,Q,nu,t,w);
if any(~far)
    h=5e-4; k=(-fix(3.5/h):fix(3.5/h));
    w=(0.5*h*pi*cosh(k*h))./cosh(0.5*pi*sinh(k*h)).^2;
    t=tanh(0.5*pi*sinh(k*h));
    keep=w>0 & abs(t)<1;
    d(~far,:)=edgeQuadrature(x2(~far),x3(~far),P,Q,nu,t(keep),w(keep));
end
end

function d=edgeQuadrature(x2,x3,P,Q,nu,t,w)
% weighted sum of the differentiated Green's function at nodes t (in
% [-1, 1]) of PQ, in the column order of edgeGradient
w=w*norm(Q-P)/2;
y2=(P(1)+Q(1))/2+t*(Q(1)-P(1))/2;
y3=(P(2)+Q(2))/2+t*(Q(2)-P(2))/2;
d=zeros(numel(x2),8);
chunk=max(1,floor(2e6/numel(t)));
col=[1 5 3 7 2 6 4 8];          % G22,2 G22,3 G32,2 G32,3 G23,2 G23,3 G33,2 G33,3
for i=1:chunk:numel(x2)
    j=i:min(numel(x2),i+chunk-1);
    g=dGreen(x2(j),x3(j),y2,y3,nu);
    for c=1:8
        d(j,c)=g{col(c)}*w';
    end
end
end

function [t,w]=gaussLegendre64()
% 64-point Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch)
persistent T W
if isempty(T)
    n=64; b=(1:n-1)./sqrt(4*(1:n-1).^2-1);
    [Vec,Dg]=eig(diag(b,1)+diag(b,-1));
    [T,i]=sort(diag(Dg)'); W=2*Vec(1,i).^2;
end
t=T; w=W;
end

function in=insideTriangle(x2,x3,A,B,C)
% points strictly inside the triangle ABC
s=@(P,Q) (Q(1)-P(1))*(x3-P(2))-(Q(2)-P(2))*(x2-P(1));
a=s(A,B); b=s(B,C); c=s(C,A);
in=(a>0&b>0&c>0)|(a<0&b<0&c<0);
end

function V=snapNearlyAlignedEdges(V)
% as in computeDisplacementPlaneStrainTriangleShearZone
tol=1e-7;
E=[1 2;2 3;3 1];
for k=1:3
    i=E(k,1); j=E(k,2);
    d=V(j,:)-V(i,:);
    L=norm(d);
    if d(2)~=0 && abs(d(2))<tol*L
        z=(V(i,2)+V(j,2))/2; V(i,2)=z; V(j,2)=z;
    end
    if d(1)~=0 && abs(d(1))<tol*L
        x=(V(i,1)+V(j,1))/2; V(i,1)=x; V(j,1)=x;
    end
end
end

% ---------------------------------------------------------------------
% Everything below is generated by derivation/generateTriangleStress.m.
% Do not edit by hand.
% ---------------------------------------------------------------------
