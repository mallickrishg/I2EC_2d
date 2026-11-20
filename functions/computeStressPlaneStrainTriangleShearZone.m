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

function [d2,d3] = dI22_v(x2,x3,m2,m3,a2,a3,t,nu)
%dI22_v
%    [D2,D3] = dI22_v(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I22_v
t2 = m2.*2.0;
t3 = m3.*2.0;
t4 = m3.*4.0;
t5 = m2.^2;
t6 = m3.^2;
t7 = nu.*4.0;
t8 = nu.^2;
t9 = t.*2.0;
t10 = t.*4.0;
t11 = t.^2;
t12 = x2.*2.0;
t13 = x3.*2.0;
t14 = x3.*4.0;
t15 = x2.^2;
t16 = x3.^2;
t19 = 1.0./pi;
t20 = m3+t+x3;
t21 = nu.*8.0;
t22 = nu-1.0;
t23 = -x2;
t25 = -x3;
t29 = m2.*x2.*-2.0;
t30 = t.*x3.*-2.0;
t31 = nu.*(3.0./2.0);
t38 = nu-1.0./2.0;
t39 = nu-3.0./4.0;
t18 = t9.*x3;
t24 = -t12;
t26 = -t13;
t27 = -t14;
t32 = t20.^2;
t33 = t7-3.0;
t35 = t20.*x3;
t36 = m2+t23;
t37 = -t31;
t40 = t9+t13;
t41 = m3+t+t25;
t42 = t21-6.0;
t43 = 1.0./t22;
t44 = t2+t24;
t45 = t36.^2;
t46 = t9+t26;
t47 = 1.0./t36;
t49 = t41.^2;
t50 = m3.*t40;
t51 = t5.*t33;
t52 = t15.*t33;
t53 = m2.*t42.*x2;
t55 = t3+t40;
t56 = t4+t10+t27;
t59 = t8+t37+5.0./8.0;
t62 = t5+t15+t29+t32;
t63 = t5+t15+t29+t35;
t48 = 1.0./t45;
t54 = m3.*t46;
t58 = t20.*t47;
t60 = t32+t45;
t64 = t45+t49;
t67 = 1.0./t62;
t69 = t36.*t63;
t70 = t6+t11+t16+t18+t45+t50;
t61 = atan(t58);
t65 = 1.0./t60;
t66 = 1.0./t64;
t68 = t67.^2;
t71 = log(t70);
t72 = t6+t11+t16+t30+t45+t54;
t73 = 1.0./t70;
t75 = 1.0./t72;
t76 = t25.*t61.*t62;
t78 = t36.*t62.*t71;
t79 = t69+t76+t78;
d2 = (t19.*t43.*(t38.*(atan(t41.*t47).*8.0-t36.*t41.*t66.*8.0)+t39.*t44.*t56.*t75))./1.6e+1-(t19.*t43.*t47.*(-t.*t33+t71.*x3+t61.*(m2.*t42+t24.*t33)+t20.*t65.*(t16.*2.0-t51-t52+t53)+t36.*t44.*t73.*x3))./8.0-(t19.*t43.*t59.*(t61.*2.0-t20.*t44.*t65+t20.*t44.*t73))./2.0-(t19.*t43.*t48.*(t61.*(t16.*2.0-t51-t52+t53)+t.*t33.*t36+t25.*t36.*t71))./8.0-(t19.*t43.*t48.*t67.*t79.*x3)./4.0+(t19.*t43.*t47.*t67.*x3.*(t63+t36.*t44+t62.*t71+t25.*t44.*t61+t36.*t44.*t71+t35.*t62.*t65+t36.*t44.*t62.*t73))./4.0-(t19.*t43.*t44.*t47.*t68.*t79.*x3)./4.0;
if nargout > 1
    d3 = (t19.*t43.*(t39.*log(t72).*4.0+t38.*t45.*t66.*8.0+t39.*t56.*t75.*(t3+t46)))./1.6e+1+(t19.*t43.*t59.*(t71+t36.*t44.*t65+t20.*t55.*t73))./2.0-(t19.*t43.*t47.*(-t36.*t71+t61.*x3.*4.0+t36.*t65.*(t16.*2.0-t51-t52+t53)+t25.*t36.*t55.*t73))./8.0-(t19.*t43.*t47.*t67.*t79)./4.0-(t19.*t43.*t47.*t67.*x3.*(-t61.*t62+t36.*(m3+t+t13)+t25.*t55.*t61+t36.*t55.*t71+t25.*t36.*t62.*t65+t36.*t55.*t62.*t73))./4.0+(t19.*t43.*t47.*t55.*t68.*t79.*x3)./4.0;
end
end

function [d2,d3] = dI22_g(x2,x3,m2,m3,a2,a3,t,nu)
%dI22_g
%    [D2,D3] = dI22_g(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I22_g
t2 = a2.*m2;
t3 = a2.*m3;
t4 = a3.*m2;
t5 = a3.*m3;
t6 = a2.*t;
t7 = a3.*t;
t8 = a2.*x2;
t9 = a2.*x3;
t10 = a3.*x2;
t11 = a3.*x3;
t12 = m3.*x3;
t13 = m3+x3;
t14 = a2.^2;
t15 = a2.^3;
t16 = a3.^2;
t18 = a3.^3;
t20 = m2.*2.0;
t21 = m2.*3.0;
t22 = m3.*2.0;
t24 = m2.^2;
t25 = m3.^2;
t26 = nu.^2;
t27 = t.^2;
t28 = x2.*2.0;
t29 = x2.*3.0;
t30 = x3.*2.0;
t31 = x2.^2;
t32 = x3.^2;
t38 = 1.0./pi;
t39 = -m2;
t40 = nu-1.0;
t41 = -x2;
t43 = -x3;
t53 = m2.*x2.*-2.0;
t61 = m2./2.0;
t62 = m2.*(3.0./4.0);
t63 = m3.*(3.0./4.0);
t64 = nu.*(3.0./2.0);
t65 = x2./2.0;
t66 = x3./2.0;
t67 = x2.*(3.0./4.0);
t77 = nu-1.0./2.0;
t78 = nu-1.0./4.0;
t79 = nu-3.0./4.0;
t80 = nu-5.0./4.0;
t17 = t14.^2;
t19 = t16.^2;
t33 = t6.*2.0;
t34 = t7.*2.0;
t36 = t12.*2.0;
t37 = t12.*3.0;
t42 = -t28;
t44 = -t29;
t45 = a2.*t13;
t46 = nu.*t13;
t47 = m3.*t22;
t48 = -t4;
t49 = t30.*x3;
t50 = -t8;
t51 = -t9;
t52 = -t11;
t55 = t3.*t15;
t56 = a2.*t6;
t57 = a3.*t7;
t58 = t11.*t18;
t59 = t13+t22;
t68 = a3.*t13.*2.0;
t70 = t7.*t11.*t16;
t71 = t6.*t13;
t72 = t7+t13;
t73 = m2+t41;
t74 = m3+t43;
t75 = -t62;
t76 = -t64;
t81 = -t65;
t84 = t25.*(3.0./4.0);
t85 = 1.0./t40;
t86 = t32.*(3.0./4.0);
t87 = t32.*(5.0./4.0);
t91 = t14+t16;
t92 = t6.^2;
t93 = t7.^2;
t103 = t13+t39+x2;
t105 = (a3.*t11)./4.0;
t106 = t14.*t40;
t123 = t16.*t77;
t138 = (t6.*t15.*t79)./2.0;
t54 = -t36;
t69 = t45.*4.0;
t82 = m3+t73;
t83 = -t46;
t89 = t13.*t34;
t90 = -t55;
t94 = t20+t42;
t95 = t21+t44;
t96 = t73.^2;
t97 = a3.*t72;
t98 = a2.*t73;
t99 = a3.*t73;
t100 = a2.*t74;
t102 = t13+t73;
t104 = t56./2.0;
t108 = a3.*t74.*2.0;
t110 = t6.*t73;
t111 = t6+t73;
t112 = t7+t74;
t113 = t34.*t79;
t115 = nu.*t103;
t121 = t33.*t73;
t124 = t27.*t91;
t126 = t13.*2.0+t34;
t128 = (m3.*t73)./3.0;
t131 = t61+t81;
t134 = t14.*t16.*t59;
t137 = t74.*t80;
t143 = t26+t76+5.0./8.0;
t144 = t3+t9+t10+t48;
t148 = t3+t10+t48+t51;
t151 = t63+t67+t75;
t160 = t106+t123;
t162 = (t15.*t73.*t79)./2.0;
t175 = t2+t5+t11+t50+t56+t57;
t181 = t2+t5+t50+t52+t56+t57;
t107 = t99.*4.0;
t109 = t13.*t83;
t114 = nu.*t102;
t116 = t96.*3.0;
t117 = -t99;
t118 = -t100;
t120 = a2.*t94;
t122 = t110.*3.0;
t125 = t15.*t99;
t127 = t97./2.0;
t129 = a3.*t112;
t130 = nu.*t96.*4.0;
t132 = t33+t94;
t136 = a2.*t18.*t95;
t139 = t3.*t14.*t111;
t140 = -t134;
t141 = t63+t83+x3;
t142 = a2.*t131;
t145 = a3.*t79.*t98;
t149 = t144.^2;
t152 = a2.*t144;
t153 = a3.*t144;
t154 = t63+t66+t83;
t161 = t148.^2;
t168 = (t16.*t78.*t98)./2.0;
t171 = t56+t97+t98;
t172 = (t18.*t79.*t112)./2.0;
t173 = t113+t137;
t174 = t82.*t151;
t178 = a2.*t16.*(t71-t128).*-3.0;
t179 = a2.*t175;
t180 = a3.*t175;
t182 = t175.^2;
t188 = t181.^2;
t189 = t32+t37+t47+t96+t110;
t209 = t25+t32+t36+t89+t92+t93+t96+t121;
t119 = -t107;
t133 = -t114;
t135 = -t130;
t146 = t45+t117;
t147 = t68+t120;
t155 = t14.*t141;
t157 = t99+t118;
t166 = t108+t120;
t176 = t12+t32+t96+t122;
t177 = t56+t98+t129;
t183 = t125.*t141.*8.0;
t185 = t36+t84+t86+t109;
t190 = t18.*t98.*t154.*8.0;
t192 = (a3.*t14.*t173)./2.0;
t196 = t104+t127+t142;
t197 = a3.*t14.*t189;
t202 = t152+t180;
t203 = t149+t182;
t208 = t161+t188;
t210 = t58+t90+t125+t136+t140;
t214 = 1.0./t209;
t150 = t69+t119;
t158 = 1.0./t146;
t164 = t.*t147;
t167 = t103.*t133;
t169 = t.*t166;
t170 = 1.0./t157;
t186 = t49+t116+t135;
t187 = t18.*t176;
t191 = -t190;
t194 = t17.*t185.*4.0;
t198 = -t197;
t200 = t105+t145+t155;
t206 = 1.0./t203;
t212 = 1.0./t208;
t215 = t214.^2;
t222 = t138+t162+t168+t172+t192;
t159 = t158.^2;
t184 = t57.*t79.*t150;
t193 = t19.*t186;
t195 = -t194;
t199 = t158.*t171;
t205 = t170.*t177;
t211 = t24+t25+t31+t32+t36+t53+t124+t164;
t217 = t24+t25+t31+t32+t53+t54+t124+t169;
t218 = t12+t87+t167+t174;
t223 = t70+t139+t178+t187+t198;
t201 = atan(t199);
t207 = atan(t205);
t213 = log(t211);
t216 = 1.0./t211;
t219 = log(t217);
t220 = 1.0./t217;
t221 = t14.*t16.*t218.*4.0;
t227 = t146.*t223;
t224 = a3.*t150.*t200.*t213;
t226 = t18.*t146.*t209.*t213;
t228 = t201.*t209.*t210;
t230 = t183+t191+t193+t195+t221;
t225 = -t224;
t229 = -t228;
t231 = t201.*t230;
t232 = t184+t225+t231;
t233 = t226+t227+t229;
et1 = (t38.*t85.*(t219.*((t15.*t79)./2.0+(a2.*t16.*t78)./2.0)+a3.*t160.*t207+t132.*t220.*t222-t157.*t160.*t212.*(a2.*t148+a3.*t181)))./2.0-(t38.*t85.*t158.*(t201.*(t19.*(m2.*-6.0+x2.*6.0+nu.*t94.*4.0)+t14.*t16.*(m2.*(3.0./2.0)+t115+t133-x2.*(3.0./2.0)).*4.0-a3.*t15.*t141.*8.0+a2.*t18.*t154.*8.0)+t7.*t16.*t79.*4.0-t16.*t200.*t213.*4.0-t202.*t206.*t230+a2.*t16.*t79.*t150.*t213+a3.*t132.*t150.*t200.*t216))./8.0-t38.*t85.*t143.*(-a3.*t201+(a2.*t213)./2.0+t132.*t196.*t216+t146.*t202.*t206)+(a3.*t38.*t85.*t159.*t232)./8.0+(t11.*t38.*t85.*t159.*t214.*t233)./4.0;
et2 = t38.*t85.*t158.*t214.*x3.*(a3.*t223-t146.*(t3.*t14+t3.*t16+t18.*(t6.*3.0+t94)-a3.*t14.*(t6+t94))+t201.*t209.*(a3.*t15+a2.*t18.*3.0)+t19.*t209.*t213+t132.*t201.*t210-t18.*t132.*t146.*t213+t202.*t206.*t209.*t210-t18.*t132.*t146.*t209.*t216).*(-1.0./4.0)-(t38.*t85.*t132.*t158.*t215.*t233.*x3)./4.0;
d2 = et1+et2;
if nargout > 1
    et3 = (t38.*t85.*(t219.*((t18.*t79)./2.0+(a3.*t14.*t80)./2.0)+t220.*t222.*(t22-t30+t34)-a2.*t160.*t207-t157.*t160.*t212.*(a3.*t148-a2.*t181)))./2.0+t38.*t85.*t143.*(a2.*t201+(a3.*t213)./2.0+t126.*t196.*t216+t146.*t206.*(t153-t179))-(t38.*t85.*t158.*(t201.*(t58.*4.0-t40.*t125.*8.0-t17.*(t22-t46.*2.0+x3.*(3.0./2.0)).*4.0+t14.*t16.*(m3-t115+t133+x3.*(5.0./2.0)).*4.0+t18.*t77.*t98.*8.0)+t6.*t16.*t79.*4.0+t206.*t230.*(t153-t179)-a2.*a3.*t200.*t213.*4.0-a3.*t150.*t213.*(t16./4.0-t106)-a3.*t126.*t150.*t200.*t216))./8.0+(a2.*t38.*t85.*t159.*t232)./8.0;
    et4 = t38.*t85.*t158.*t214.*t233.*(-1.0./4.0)+(t9.*t38.*t85.*t159.*t214.*t233)./4.0-(t38.*t85.*t158.*t214.*x3.*(a2.*t223+t146.*(t18.*(t13+x3)+t7.*t18-t16.*t56.*3.0-a3.*t14.*(m3+t13.*2.0))-t201.*t209.*(t19-t14.*t16)-t126.*t201.*t210+a2.*t18.*t209.*t213+t18.*t126.*t146.*t213-t206.*t209.*t210.*(t153-t179)+t18.*t126.*t146.*t209.*t216))./4.0+(t38.*t85.*t126.*t158.*t215.*t233.*x3)./4.0;
    d3 = et3+et4;
end
end

function [d2,d3] = dI33_v(x2,x3,m2,m3,a2,a3,t,nu)
%dI33_v
%    [D2,D3] = dI33_v(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I33_v
t2 = m3+t;
t3 = m2.*2.0;
t6 = m2.^2;
t7 = m3.^2;
t8 = nu.*2.0;
t9 = nu.^2;
t10 = t.*2.0;
t12 = t.^2;
t13 = x2.*2.0;
t14 = x3.*2.0;
t15 = x3.*4.0;
t16 = x2.^2;
t17 = x3.^2;
t20 = 1.0./pi;
t22 = nu.*8.0;
t23 = nu-1.0;
t24 = -x2;
t26 = -x3;
t31 = m2.*x2.*-2.0;
t32 = t.*x3.*-2.0;
t34 = nu.*(3.0./2.0);
t39 = nu-3.0./4.0;
t19 = t10.*x3;
t21 = t2+x3;
t25 = -t13;
t27 = -t14;
t28 = -t15;
t29 = t2.*x3;
t33 = t2.^2;
t37 = m2+t24;
t38 = -t34;
t40 = t10+t14;
t41 = t2+t26;
t42 = 1.0./t23;
t43 = t2.*2.0;
t47 = t8-3.0./2.0;
t48 = t22-8.0;
t53 = t2.*t14;
t56 = t16.*t39.*4.0;
t58 = t6.*t39.*4.0;
t35 = t21.^2;
t44 = t3+t25;
t45 = t37.^2;
t46 = t10+t27;
t49 = 1.0./t37;
t51 = t41.^2;
t52 = m3.*t40;
t55 = t14+t43;
t57 = t2.*4.0+t28;
t59 = m2.*t47.*x2.*4.0;
t62 = t9+t38+5.0./8.0;
t70 = t6+t16+t17+t31+t33+t53;
t50 = 1.0./t45;
t54 = m3.*t46;
t61 = t21.*t49;
t63 = t35+t45;
t65 = t17+t29+t45;
t66 = t45+t51;
t71 = t7+t12+t17+t19+t45+t52;
t73 = 1.0./t70;
t64 = atan(t61);
t67 = 1.0./t63;
t68 = 1.0./t66;
t69 = t37.*t65;
t72 = log(t71);
t74 = t73.^2;
t75 = t7+t12+t17+t32+t45+t54;
t76 = 1.0./t71;
t77 = 1.0./t75;
t79 = t64.*t70.*x3;
t80 = t69+t79;
d2 = (t20.*t42.*(t48.*atan(t41.*t49)-t37.*t41.*t48.*t68+t39.*t44.*t57.*t77))./1.6e+1-(t20.*t42.*t62.*(t64.*2.0-t21.*t44.*t67+t21.*t44.*t76))./2.0+(t20.*t42.*t49.*(t72.*x3+t64.*(m2.*t47.*4.0-t39.*x2.*8.0)+t21.*t67.*(t17.*2.0-t56-t58+t59)+t37.*t44.*t76.*x3))./8.0+(t20.*t42.*t50.*(t64.*(t17.*2.0-t56-t58+t59)+t26.*t37.*t72))./8.0-(t20.*t42.*t50.*t73.*t80.*x3)./4.0+(t20.*t42.*t49.*t73.*x3.*(t65+t37.*t44+t44.*t64.*x3+t21.*t26.*t67.*t70))./4.0-(t20.*t42.*t44.*t49.*t74.*t80.*x3)./4.0;
if nargout > 1
    d3 = (t20.*t42.*(t39.*log(t75).*4.0+t45.*t48.*t68+t39.*t57.*t77.*(t27+t43)))./1.6e+1+(t20.*t42.*t62.*(t72+t37.*t44.*t67+t21.*t55.*t76))./2.0+(t20.*t42.*t49.*(-t37.*t72+t64.*x3.*4.0+t37.*t67.*(t17.*2.0-t56-t58+t59)+t26.*t37.*t55.*t76))./8.0-(t20.*t42.*t49.*t73.*t80)./4.0-(t20.*t42.*t49.*t73.*x3.*(t37.*(t2+t14)+t64.*t70+t55.*t64.*x3+t37.*t67.*t70.*x3))./4.0+(t20.*t42.*t49.*t55.*t74.*t80.*x3)./4.0;
end
end

function [d2,d3] = dI33_g(x2,x3,m2,m3,a2,a3,t,nu)
%dI33_g
%    [D2,D3] = dI33_g(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I33_g
t2 = a2.*m2;
t3 = a2.*m3;
t4 = a3.*m2;
t5 = a3.*m3;
t6 = a2.*t;
t7 = a3.*t;
t8 = a2.*x2;
t9 = a2.*x3;
t10 = a3.*x2;
t11 = a3.*x3;
t12 = m3.*x3;
t13 = m3+x3;
t14 = a2.^2;
t15 = a2.^3;
t16 = a3.^2;
t18 = a3.^3;
t20 = m2.*2.0;
t21 = m3.*2.0;
t22 = m2.*6.0;
t23 = m2.^2;
t24 = m3.^2;
t25 = nu.*4.0;
t26 = nu.^2;
t27 = t.^2;
t28 = x2.*2.0;
t29 = x3.*2.0;
t30 = x2.*6.0;
t31 = x2.^2;
t32 = x3.^2;
t39 = 1.0./pi;
t40 = -m2;
t41 = nu-1.0;
t42 = -x2;
t44 = -x3;
t56 = m2.*x2.*-2.0;
t62 = m2./2.0;
t63 = m3./3.0;
t64 = m2.*(3.0./4.0);
t65 = m3.*(3.0./4.0);
t66 = nu.*(3.0./2.0);
t67 = x2./2.0;
t68 = x3./2.0;
t69 = x2.*(3.0./4.0);
t78 = nu-1.0./2.0;
t79 = nu-1.0./4.0;
t80 = nu-3.0./4.0;
t81 = nu-5.0./4.0;
t17 = t14.^2;
t19 = t16.^2;
t33 = t6.*2.0;
t34 = t7.*2.0;
t36 = t12.*2.0;
t37 = t12.*3.0;
t38 = t12.*4.0;
t43 = -t28;
t45 = a2.*t13;
t46 = nu.*t13;
t47 = m3.*t21;
t48 = t24.*3.0;
t49 = -t4;
t50 = t29.*x3;
t51 = t32.*3.0;
t52 = t32.*5.0;
t53 = -t8;
t54 = -t9;
t55 = -t11;
t58 = t12.*8.0;
t59 = a2.*t6;
t60 = a3.*t7;
t70 = a3.*t13.*2.0;
t71 = t7.*t11.*t16;
t72 = t25-3.0;
t73 = t7+t13;
t74 = m2+t42;
t75 = m3+t44;
t76 = -t64;
t77 = -t66;
t82 = -t67;
t85 = 1.0./t41;
t86 = t6.*t13.*3.0;
t88 = t63+x3;
t89 = t14+t16;
t90 = t6.^2;
t91 = t7.^2;
t102 = t13+t40+x2;
t103 = (t3.*t15)./3.0;
t105 = (a3.*t11)./4.0;
t106 = (t11.*t18)./3.0;
t107 = t16.*t41;
t121 = t14.*t78;
t136 = (t6.*t15.*t80)./2.0;
t57 = -t36;
t83 = m3+t74;
t84 = -t46;
t87 = t13.*t34;
t92 = t20+t43;
t93 = t74.^2;
t94 = a3.*t73;
t95 = a2.*t74;
t96 = a3.*t74;
t97 = a2.*t75;
t98 = m3.*t74;
t101 = t13+t74;
t104 = t59./2.0;
t108 = a3.*t75.*2.0;
t109 = t13.*t46.*-4.0;
t110 = -t103;
t111 = t6.*t74;
t112 = t6+t74;
t113 = t7+t75;
t114 = t34.*t80;
t119 = t33.*t74;
t122 = t27.*t89;
t123 = t13.*2.0+t34;
t126 = t25.*t102;
t129 = t62+t82;
t133 = t75.*t79;
t138 = t14.*t16.*t88;
t141 = t26+t77+5.0./8.0;
t143 = t3+t9+t10+t49;
t147 = t3+t10+t49+t54;
t150 = t65+t69+t76;
t158 = t107+t121;
t160 = (t15.*t74.*t80)./2.0;
t173 = t2+t5+t11+t53+t59+t60;
t180 = t2+t5+t53+t55+t59+t60;
t115 = t93.*3.0;
t116 = -t96;
t117 = -t97;
t118 = a2.*t92;
t120 = t111.*3.0;
t124 = t94./2.0;
t127 = a3.*t113;
t130 = t33+t92;
t131 = nu.*t101.*-4.0;
t132 = nu.*t93.*-4.0;
t134 = (t18.*t95)./3.0;
t135 = (t15.*t96)./3.0;
t137 = t3.*t14.*t112;
t139 = t65+t84+x3;
t140 = a2.*t129;
t144 = a3.*t80.*t95;
t149 = t143.^2;
t151 = a2.*t143;
t152 = a3.*t143;
t153 = t65+t68+t84;
t159 = t147.^2;
t166 = (t16.*t81.*t95)./2.0;
t169 = t59+t94+t95;
t170 = (t18.*t80.*t113)./2.0;
t171 = -a2.*t16.*(t86-t98);
t172 = t114+t133;
t175 = t48+t51+t58+t109;
t177 = t83.*t150.*4.0;
t178 = a2.*t173;
t179 = a3.*t173;
t182 = t173.^2;
t187 = t180.^2;
t188 = t32+t37+t47+t93+t111;
t210 = t24+t32+t36+t87+t90+t91+t93+t119;
t142 = -t134;
t145 = t45+t116;
t146 = t70+t118;
t154 = t14.*t139;
t155 = t96+t117;
t164 = t108+t118;
t165 = t102.*t131;
t174 = t12+t32+t93+t120;
t176 = t59+t95+t127;
t183 = t15.*t96.*t139.*8.0;
t184 = t50+t115+t132;
t186 = t17.*t175;
t189 = t18.*t95.*t153.*8.0;
t192 = (a3.*t14.*t172)./2.0;
t194 = t104+t124+t140;
t195 = a3.*t14.*t188;
t201 = t151+t179;
t203 = t149+t182;
t209 = t159+t187;
t214 = 1.0./t210;
t156 = 1.0./t145;
t163 = t.*t146;
t167 = t.*t164;
t168 = 1.0./t155;
t181 = t59.*t72.*t145;
t185 = t18.*t174;
t190 = -t186;
t191 = -t189;
t193 = t19.*t184;
t196 = -t195;
t198 = t105+t144+t154;
t207 = 1.0./t203;
t212 = 1.0./t209;
t215 = t214.^2;
t218 = t38+t52+t165+t177;
t221 = t106+t110+t135+t138+t142;
t223 = t136+t160+t166+t170+t192;
t157 = t156.^2;
t197 = t156.*t169;
t200 = t23+t24+t27+t31+t32+t36+t56+t163;
t206 = t168.*t176;
t211 = t23+t24+t31+t32+t36+t56+t122+t163;
t217 = t23+t24+t31+t32+t56+t57+t122+t167;
t222 = t14.*t16.*t218;
t224 = t71+t137+t171+t185+t196;
t199 = atan(t197);
t202 = log(t200);
t205 = 1.0./t200;
t208 = atan(t206);
t213 = log(t211);
t216 = 1.0./t211;
t219 = log(t217);
t220 = 1.0./t217;
t228 = (t145.*t224)./4.0;
t230 = t183+t190+t191+t193+t222;
t225 = a3.*t145.*t198.*t213.*4.0;
t226 = (a3.*t14.*t145.*t210.*t213)./4.0;
t229 = t199.*t210.*t221.*(3.0./4.0);
t231 = t199.*t230;
t227 = -t226;
t232 = -t231;
t233 = t181+t225+t232;
t234 = t227+t228+t229;
et1 = (t39.*t85.*(t219.*((t15.*t80)./2.0+(a2.*t16.*t81)./2.0)+a3.*t158.*t208+t130.*t220.*t223-t155.*t158.*t212.*(a2.*t147+a3.*t180)))./2.0+(t39.*t85.*t156.*(t199.*(t19.*(-t22+t30+t25.*t92)-a3.*t15.*t139.*8.0+a2.*t18.*t153.*8.0+t14.*t16.*(t22-t30+t126+t131))-a3.*t59.*t72-t16.*t198.*t213.*4.0-t201.*t207.*t230+a2.*t16.*t80.*t145.*t213.*4.0+a3.*t130.*t145.*t198.*t216.*4.0))./8.0-t39.*t85.*t141.*(-a3.*t199+(a2.*t202)./2.0+t130.*t194.*t205+t145.*t201.*t207)+(a3.*t39.*t85.*t157.*t233)./8.0+t11.*t39.*t85.*t157.*t214.*t234;
et2 = -t39.*t44.*t85.*t156.*t214.*(a3.*t224.*(-1.0./4.0)+(t145.*(t3.*t14+t3.*t16+t18.*(t6.*3.0+t92)-a3.*t14.*(t6+t92)))./4.0+t199.*t210.*((a3.*t15)./3.0-(a2.*t18)./3.0).*(3.0./4.0)+t130.*t199.*t221.*(3.0./4.0)+(t14.*t16.*t210.*t213)./4.0+t201.*t207.*t210.*t221.*(3.0./4.0)-(a3.*t14.*t130.*t145.*t213)./4.0-(a3.*t14.*t130.*t145.*t210.*t216)./4.0)-t39.*t85.*t130.*t156.*t215.*t234.*x3;
d2 = et1+et2;
if nargout > 1
    et3 = (t39.*t85.*(t219.*((t18.*t80)./2.0+(a3.*t14.*t79)./2.0)+t220.*t223.*(t21-t29+t34)-a2.*t158.*t208-t155.*t158.*t212.*(a3.*t147-a2.*t180)))./2.0+t39.*t85.*t141.*(a2.*t199+(a3.*t202)./2.0+t123.*t194.*t205+t145.*t207.*(t152-t178))-(t39.*t85.*t156.*(-t199.*(t11.*t18.*4.0-t17.*(t13.*6.0+t21-t46.*8.0)+t14.*t16.*(t13.*4.0+t131+x3.*6.0-nu.*t102.*4.0)-t15.*t41.*t96.*8.0+t18.*t78.*t95.*8.0)+t6.*t14.*t72-t207.*t230.*(t152-t178)+a2.*a3.*t198.*t213.*4.0+a3.*t145.*t213.*(t16./4.0-t14.*t41).*4.0+a3.*t123.*t145.*t198.*t216.*4.0))./8.0+(a2.*t39.*t85.*t157.*t233)./8.0;
    et4 = -t39.*t85.*t156.*t214.*t234+t9.*t39.*t85.*t157.*t214.*t234-t39.*t85.*t156.*t214.*x3.*((a2.*t224)./4.0+(t145.*(t18.*(t13+x3)+t7.*t18-t16.*t59.*3.0-a3.*t14.*(m3+t13.*2.0)))./4.0+t123.*t199.*t221.*(3.0./4.0)+t199.*t210.*(t19./3.0+t14.*t16).*(3.0./4.0)-(a3.*t15.*t210.*t213)./4.0+t207.*t210.*t221.*(t152-t178).*(3.0./4.0)-(a3.*t14.*t123.*t145.*t213)./4.0-(a3.*t14.*t123.*t145.*t210.*t216)./4.0)-t39.*t44.*t85.*t123.*t156.*t215.*t234;
    d3 = et3+et4;
end
end

function [d2,d3] = dI23_v(x2,x3,m2,m3,a2,a3,t,nu)
%dI23_v
%    [D2,D3] = dI23_v(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I23_v
t2 = m3+t;
t3 = m2.*2.0;
t5 = m3.^2;
t6 = t.*2.0;
t7 = t.^2;
t8 = x2.*2.0;
t9 = x3.*2.0;
t10 = x3.^2;
t12 = 1.0./pi;
t14 = nu-1.0;
t15 = -x2;
t18 = t.*x3.*-2.0;
t19 = m2./2.0;
t20 = x2./2.0;
t23 = nu-1.0./2.0;
t24 = nu-3.0./4.0;
t11 = t6.*x3;
t13 = t2+x3;
t16 = -t8;
t17 = -t9;
t22 = m2+t15;
t25 = -t20;
t26 = t6+t9;
t29 = 1.0./t14;
t37 = t2.*2.0+t9;
t21 = t13.^2;
t27 = 1.0./t13;
t30 = t3+t16;
t31 = t22.^2;
t32 = t6+t17;
t33 = 1.0./t22;
t34 = m3.*t26;
t35 = t2.*t22;
t38 = t19+t25;
t28 = 1.0./t21;
t36 = m3.*t32;
t39 = t13.*t33;
t40 = t22.*t27;
t41 = t21+t31;
t44 = t5+t7+t10+t11+t31+t34;
t42 = atan(t39);
t43 = 1.0./t41;
t45 = t5+t7+t10+t18+t31+t36;
t46 = 1.0./t44;
t47 = t46.^2;
t48 = 1.0./t45;
t49 = t42.*t44;
t50 = -t49;
t51 = t35+t50;
d2 = (t12.*t29.*log(t45))./1.6e+1-t12.*t23.*(-log(-t40)+log(t28.*t44)./2.0+t21.*t43+t30.*t38.*t46-1.0)-(t12.*t24.*t29.*(log(t44)+t22.*t30.*t46+t13.*t43.*x3.*4.0))./4.0-(t12.*t29.*t46.*x3.*(t2-t30.*t42+t13.*t43.*t44))./4.0+(t12.*t22.*t29.*t30.*t48)./1.6e+1+(t12.*t29.*t30.*t47.*t51.*x3)./4.0;
if nargout > 1
    d3 = t12.*t23.*(t40+atan(t40)-t13.*t22.*t43+t21.*t38.*t46.*(t28.*t37-t27.^3.*t44.*2.0))+(t12.*t29.*t46.*t51)./4.0-(t12.*t24.*t29.*(t42.*4.0-t22.*t37.*t46+t22.*t43.*x3.*4.0))./4.0-(t12.*t29.*t46.*x3.*(t37.*t42+t22.*t43.*t44))./4.0+(t12.*t22.*t29.*t48.*(t2.*2.0+t17))./1.6e+1-(t12.*t29.*t37.*t47.*t51.*x3)./4.0;
end
end

function [d2,d3] = dI23_g(x2,x3,m2,m3,a2,a3,t,nu)
%dI23_g
%    [D2,D3] = dI23_g(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I23_g
t2 = a2+a3;
t3 = a2.*m2;
t4 = a2.*m3;
t5 = a3.*m2;
t6 = a3.*m3;
t7 = a2.*t;
t8 = a3.*t;
t9 = a2.*x2;
t10 = a2.*x3;
t11 = a3.*x2;
t12 = a3.*x3;
t13 = m3+x3;
t14 = a2.^2;
t15 = a2.^3;
t16 = a3.^2;
t17 = a3.^3;
t19 = m2.*2.0;
t20 = m3.*2.0;
t22 = m2.^2;
t23 = m3.^2;
t24 = t.^2;
t25 = x2.*2.0;
t26 = x3.*2.0;
t27 = x2.^2;
t28 = x3.^2;
t33 = 1.0./pi;
t34 = -a3;
t35 = 1.0./a2;
t36 = 1.0./a3;
t37 = nu-1.0;
t38 = -x2;
t40 = -x3;
t46 = m2.*x2.*-2.0;
t47 = m3.*x3.*-2.0;
t53 = m2./4.0;
t54 = m3./3.0;
t55 = m3./4.0;
t56 = x2./4.0;
t57 = x3.*(3.0./4.0);
t59 = m3.*x3.*(3.0./2.0);
t66 = nu-1.0./2.0;
t67 = nu-3.0./4.0;
t29 = t7.*2.0;
t30 = t8.*2.0;
t32 = t20.*x3;
t39 = -t25;
t41 = a2.*t13;
t42 = -t5;
t43 = -t9;
t44 = -t10;
t45 = -t12;
t48 = a2.*t4;
t49 = a2.*t7;
t50 = a3.*t8;
t51 = a3.*t12;
t52 = t13+t20;
t58 = a3.*t13.*2.0;
t60 = t8+t13;
t61 = a2+t34;
t62 = a2.*t36;
t63 = a3.*t35;
t64 = m2+t38;
t65 = m3+t40;
t68 = -t56;
t69 = t23./2.0;
t70 = 1.0./t37;
t73 = t54+x3;
t74 = t14+t16;
t76 = t7.^2;
t77 = t8.^2;
t85 = t4.*t13.*t14;
t100 = t55+t57;
t71 = t13.*t30;
t72 = -t48;
t75 = atan(t63);
t78 = t19+t39;
t79 = t64.^2;
t80 = a3.*t60;
t81 = a2.*t64;
t82 = a3.*t64;
t83 = a2.*t65;
t84 = m3.*t64;
t86 = t7.*t52;
t87 = t60.^2;
t88 = a3.*t65.*2.0;
t90 = t7+t64;
t91 = t8+t65;
t92 = 1.0./t60;
t94 = 1.0./t74;
t96 = t34.*t64;
t99 = t29.*t64;
t101 = t24.*t74;
t102 = t13.*2.0+t30;
t104 = t8.*t17.*t64;
t106 = t7.*t73.*3.0;
t107 = t53+t68;
t109 = t7.*t64.*(3.0./2.0);
t111 = (t15.*t65)./4.0;
t114 = t4+t10+t11+t42;
t116 = a2.*t16.*t100;
t119 = t4+t11+t42+t44;
t139 = t3+t6+t12+t43+t49+t50;
t145 = t3+t6+t43+t45+t49+t50;
t89 = a3.*t81;
t93 = 1.0./t87;
t95 = t94.^2;
t97 = -t83;
t98 = a2.*t78;
t103 = t90.^2;
t105 = a3.*t91;
t108 = t29+t78;
t110 = -t106;
t112 = (t14.*t82)./4.0;
t115 = t17.*t107;
t117 = t41+t96;
t120 = t114.^2;
t121 = t84+t86;
t122 = a2.*t114;
t123 = a3.*t114;
t124 = -t116;
t127 = t119.^2;
t137 = t49+t80+t81;
t143 = a2.*t139;
t144 = a3.*t139;
t147 = t139.^2;
t148 = t36.*t92.*t114;
t149 = t145.^2;
t151 = t28+t59+t69+t79+t109;
t169 = t23+t28+t32+t71+t76+t77+t79+t99;
t113 = -t112;
t118 = t58+t98;
t125 = t82+t97;
t126 = 1.0./t117;
t130 = t87+t103;
t131 = t88+t98;
t132 = t84+t110;
t133 = a3.*t14.*t121;
t134 = t51+t72+t89;
t141 = t49+t81+t105;
t142 = t92.*t117;
t150 = -t148;
t156 = a2.*t16.*t151.*2.0;
t161 = t122+t144;
t162 = t120+t147;
t168 = t127+t149;
t173 = 1.0./t169;
t129 = t.*t118;
t135 = t.*t131;
t136 = 1.0./t125;
t138 = 1.0./t130;
t140 = t17.*t132;
t146 = log(t142);
t152 = t62+t150;
t154 = t126.*t137;
t159 = -t156;
t166 = 1.0./t162;
t171 = 1.0./t168;
t174 = t173.^2;
t179 = t111+t113+t115+t124;
t153 = t152.^2;
t158 = atan(t154);
t165 = t136.*t141;
t170 = t22+t23+t27+t28+t32+t46+t101+t129;
t176 = t22+t23+t27+t28+t46+t47+t101+t135;
t155 = t153+1.0;
t163 = a2.*a3.*t158;
t167 = atan(t165);
t172 = log(t170);
t175 = 1.0./t170;
t177 = log(t176);
t178 = 1.0./t176;
t181 = t2.*t34.*t61.*t158.*t169;
t157 = log(t155);
t160 = 1.0./t155;
t182 = a2.*t16.*t169.*t172;
t183 = -t182;
t184 = t85+t104+t133+t140+t159+t181+t183;
d2 = t33.*t70.*(a2.*a3.*(a3.*t167-t125.*t171.*(a2.*t119+a3.*t145))+(a3.*t2.*t61.*t177)./4.0+(t2.*t61.*t108.*t125.*t178)./4.0).*(-1.0./4.0)-t33.*t67.*t70.*(t172.*(t17./4.0-(a3.*t14)./4.0)+a3.*(t163+t134.*t161.*t166)+t108.*t175.*t179)+t33.*t36.*t66.*t94.*(t16+t16.*t146-(t16.*t157)./2.0+a2.*a3.*t75+a2.*t34.*t158-t74.*t87.*t138+a2.*t117.*t161.*t166+a3.*t142.*t152.*t160)+(t33.*t70.*t95.*t173.*x3.*(-t6.*t16-t8.*t17+t34.*t48+a2.*t16.*(t7.*(3.0./2.0)+t78).*2.0+a2.*t16.*t108.*t172+a3.*t2.*t61.*t108.*t158+a2.*t16.*t108.*t169.*t175+a3.*t2.*t61.*t161.*t166.*t169))./4.0+(t33.*t70.*t95.*t108.*t174.*t184.*x3)./4.0;
if nargout > 1
    et1 = (t33.*t70.*(a2.*a3.*(a2.*t167-t125.*t171.*(a2.*t145+t34.*t119))+(a2.*t2.*t61.*t177)./4.0-(t2.*t61.*t125.*t178.*(t20-t26+t30))./4.0))./4.0+t33.*t67.*t70.*(-t172.*(t15./4.0+a2.*t16.*(3.0./4.0))+a3.*(t16.*t158+t134.*t166.*(t123-t143))+t102.*t175.*t179)+t33.*t36.*t66.*t94.*(t14.*t75-t14.*t158+t74.*atan(t90.*t92)+t80.*(a2.*t92-t93.*t117)+a2.*a3.*t146-(a2.*a3.*t157)./2.0-t60.*t74.*t90.*t138-a2.*t117.*t166.*(t123-t143)+a3.*t117.*t152.*t160.*(t62.*t92-t36.*t93.*t114))+(t33.*t70.*t95.*t173.*t184)./4.0;
    et2 = t33.*t70.*t95.*t173.*x3.*(-t4.*t14+t7.*t17.*3.0+t7.*t14.*t34+a2.*t16.*(m3.*(3.0./2.0)+t26).*2.0+a2.*t16.*t102.*t172+a3.*t2.*t61.*t102.*t158+a2.*t16.*t102.*t169.*t175+a3.*t2.*t61.*t166.*t169.*(t123-t143)).*(-1.0./4.0)-(t33.*t70.*t95.*t102.*t174.*t184.*x3)./4.0;
    d3 = et1+et2;
end
end

function [d2,d3] = dI23_h(x2,x3,m2,m3,a2,a3,t,nu)
%dI23_h
%    [D2,D3] = dI23_h(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I23_h
t2 = a2+a3;
t3 = a2.*m2;
t4 = a2.*m3;
t5 = a3.*m2;
t6 = a3.*m3;
t7 = a3.*t;
t8 = a2.*x2;
t9 = a2.*x3;
t10 = a3.*x2;
t11 = a3.*x3;
t12 = m3+x3;
t13 = a2.^2;
t14 = a2.^3;
t15 = a3.^2;
t16 = a3.^3;
t18 = m2.*2.0;
t19 = m3.*2.0;
t21 = m2.^2;
t22 = m3.^2;
t23 = t.^2;
t24 = x2.*2.0;
t25 = x3.*2.0;
t26 = x2.^2;
t27 = x3.^2;
t28 = a2.*t.*2.0;
t32 = 1.0./pi;
t33 = -a3;
t34 = nu-1.0;
t35 = -x2;
t37 = -x3;
t43 = m2.*x2.*-2.0;
t44 = m3.*x3.*-2.0;
t51 = m2./4.0;
t52 = m3./3.0;
t53 = m3./4.0;
t54 = x2./4.0;
t55 = x3.*(3.0./4.0);
t57 = m3.*x3.*(3.0./2.0);
t63 = nu-1.0./2.0;
t64 = nu-3.0./4.0;
t29 = t7.*2.0;
t31 = t19.*x3;
t36 = -t24;
t38 = a2.*t12;
t39 = -t5;
t40 = -t8;
t41 = -t9;
t42 = -t11;
t45 = a2.*t4;
t46 = t.*t13;
t47 = a3.*t7;
t48 = a3.*t11;
t49 = t12+t19;
t50 = t12.^2;
t56 = a3.*t12.*2.0;
t58 = t7+t12;
t59 = a2+t33;
t60 = m2+t35;
t61 = m3+t37;
t65 = -t54;
t67 = t22./2.0;
t68 = 1.0./t34;
t71 = t52+x3;
t72 = t13+t15;
t73 = t13.*t23;
t74 = t7.^2;
t82 = t4.*t12.*t13;
t93 = t53+t55;
t62 = 1.0./t50;
t66 = t+t60;
t69 = t12.*t29;
t70 = -t45;
t75 = t18+t36;
t76 = t60.^2;
t77 = a3.*t58;
t78 = a2.*t60;
t79 = a3.*t60;
t80 = a2.*t61;
t81 = m3.*t60;
t83 = a2.*t.*t49;
t84 = a3.*t61.*2.0;
t87 = t7+t61;
t88 = 1.0./t72.^2;
t89 = t33.*t60;
t92 = t28.*t60;
t94 = t23.*t72;
t95 = t12.*2.0+t29;
t96 = t7.*t16.*t60;
t98 = a2.*t.*t71.*3.0;
t99 = t51+t65;
t103 = (t14.*t61)./4.0;
t106 = t4+t9+t10+t39;
t109 = a2.*t15.*t93;
t112 = t4+t10+t39+t41;
t133 = t3+t6+t11+t40+t46+t47;
t138 = t3+t6+t40+t42+t46+t47;
t85 = t66.^2;
t86 = a3.*t78;
t90 = -t80;
t91 = a2.*t75;
t97 = a3.*t87;
t100 = t28+t75;
t101 = t.*t78.*(3.0./2.0);
t102 = -t98;
t104 = (t13.*t79)./4.0;
t108 = t16.*t99;
t110 = t38+t89;
t113 = t106.^2;
t114 = t81+t83;
t115 = a2.*t106;
t116 = a3.*t106;
t117 = -t109;
t121 = t112.^2;
t131 = t46+t77+t78;
t136 = a2.*t133;
t137 = a3.*t133;
t139 = t133.^2;
t140 = t138.^2;
t153 = t22+t27+t31+t69+t73+t74+t76+t92;
t105 = -t104;
t107 = t50+t85;
t111 = t56+t91;
t119 = t79+t90;
t120 = 1.0./t110;
t124 = t84+t91;
t125 = t22+t27+t31+t85;
t126 = t81+t102;
t127 = a3.*t13.*t114;
t128 = t48+t70+t86;
t135 = t46+t78+t97;
t141 = t27+t57+t67+t76+t101;
t146 = t115+t137;
t147 = t113+t139;
t152 = t121+t140;
t157 = 1.0./t153;
t118 = 1.0./t107;
t123 = t.*t111;
t129 = t.*t124;
t130 = 1.0./t119;
t132 = 1.0./t125;
t134 = t16.*t126;
t142 = t120.*t131;
t143 = a2.*t15.*t141.*2.0;
t150 = 1.0./t147;
t155 = 1.0./t152;
t158 = t157.^2;
t163 = t103+t105+t108+t117;
t144 = atan(t142);
t145 = -t143;
t149 = t130.*t135;
t154 = t21+t22+t26+t27+t31+t43+t94+t123;
t160 = t21+t22+t26+t27+t43+t44+t94+t129;
t151 = atan(t149);
t156 = log(t154);
t159 = 1.0./t154;
t161 = log(t160);
t162 = 1.0./t160;
t165 = t2.*t33.*t59.*t144.*t153;
t166 = a2.*t15.*t153.*t156;
t167 = -t166;
t168 = t82+t96+t127+t134+t145+t165+t167;
d2 = t32.*t68.*(a2.*a3.*(a3.*t151-t119.*t155.*(a2.*t112+a3.*t138))+(a3.*t2.*t59.*t161)./4.0+(t2.*t59.*t100.*t119.*t162)./4.0).*(-1.0./4.0)-t32.*t63.*(atan(t66./t12)+t12.*t66.*t118-(t12.*t132.*(t.*2.0+t75))./2.0)-t32.*t64.*t68.*(t156.*(t16./4.0-(a3.*t13)./4.0)+a3.*(a2.*a3.*t144+t128.*t146.*t150)+t100.*t159.*t163)+(t32.*t68.*t88.*t157.*x3.*(-t6.*t15-t7.*t16+t33.*t45+a2.*t15.*(t75+a2.*t.*(3.0./2.0)).*2.0+a2.*t15.*t100.*t156+a3.*t2.*t59.*t100.*t144+a2.*t15.*t100.*t153.*t159+a3.*t2.*t59.*t146.*t150.*t153))./4.0+(t32.*t68.*t88.*t100.*t158.*t168.*x3)./4.0;
if nargout > 1
    et1 = -t32.*t63.*(log(t62.*t125)./2.0+t85.*t118-(t12.^3.*t132.*(1.0./t12.^3.*t125.*2.0-2.0./t12))./2.0)+(t32.*t68.*(a2.*a3.*(a2.*t151-t119.*t155.*(a2.*t138+t33.*t112))+(a2.*t2.*t59.*t161)./4.0-(t2.*t59.*t119.*t162.*(t19-t25+t29))./4.0))./4.0+t32.*t64.*t68.*(-t156.*(t14./4.0+a2.*t15.*(3.0./4.0))+a3.*(t15.*t144+t128.*t150.*(t116-t136))+t95.*t159.*t163)+(t32.*t68.*t88.*t157.*t168)./4.0;
    et2 = t32.*t68.*t88.*t157.*x3.*(-t4.*t13-t7.*t14+a2.*t7.*t15.*3.0+a2.*t15.*(m3.*(3.0./2.0)+t25).*2.0+a2.*t15.*t95.*t156+a3.*t2.*t59.*t95.*t144+a2.*t15.*t95.*t153.*t159+a3.*t2.*t59.*t150.*t153.*(t116-t136)).*(-1.0./4.0)-(t32.*t68.*t88.*t95.*t158.*t168.*x3)./4.0;
    d3 = et1+et2;
end
end

function [d2,d3] = dI32_v(x2,x3,m2,m3,a2,a3,t,nu)
%dI32_v
%    [D2,D3] = dI32_v(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I32_v
t2 = m3+t;
t3 = m2.*2.0;
t5 = m3.^2;
t6 = t.*2.0;
t7 = t.^2;
t8 = x2.*2.0;
t9 = x3.*2.0;
t10 = x3.^2;
t12 = 1.0./pi;
t14 = nu-1.0;
t15 = -x2;
t18 = t.*x3.*-2.0;
t19 = m2./2.0;
t20 = x2./2.0;
t23 = nu-1.0./2.0;
t24 = nu-3.0./4.0;
t11 = t6.*x3;
t13 = t2+x3;
t16 = -t8;
t17 = -t9;
t22 = m2+t15;
t25 = -t20;
t26 = t6+t9;
t29 = 1.0./t14;
t37 = t2.*2.0+t9;
t21 = t13.^2;
t27 = 1.0./t13;
t30 = t3+t16;
t31 = t22.^2;
t32 = t6+t17;
t33 = 1.0./t22;
t34 = m3.*t26;
t35 = t2.*t22;
t38 = t19+t25;
t28 = 1.0./t21;
t36 = m3.*t32;
t39 = t13.*t33;
t40 = t22.*t27;
t41 = t21+t31;
t44 = t5+t7+t10+t11+t31+t34;
t42 = atan(t39);
t43 = 1.0./t41;
t45 = t5+t7+t10+t18+t31+t36;
t46 = 1.0./t44;
t47 = t46.^2;
t48 = 1.0./t45;
t49 = t42.*t44;
t50 = -t49;
t51 = t35+t50;
d2 = (t12.*t29.*log(t45))./1.6e+1+t12.*t23.*(-log(-t40)+log(t28.*t44)./2.0+t21.*t43+t30.*t38.*t46-1.0)-(t12.*t24.*t29.*(log(t44)+t22.*t30.*t46+t13.*t43.*x3.*4.0))./4.0+(t12.*t29.*t46.*x3.*(t2-t30.*t42+t13.*t43.*t44))./4.0+(t12.*t22.*t29.*t30.*t48)./1.6e+1-(t12.*t29.*t30.*t47.*t51.*x3)./4.0;
if nargout > 1
    d3 = -t12.*t23.*(t40+atan(t40)-t13.*t22.*t43+t21.*t38.*t46.*(t28.*t37-t27.^3.*t44.*2.0))-(t12.*t29.*t46.*t51)./4.0-(t12.*t24.*t29.*(t42.*4.0-t22.*t37.*t46+t22.*t43.*x3.*4.0))./4.0+(t12.*t29.*t46.*x3.*(t37.*t42+t22.*t43.*t44))./4.0+(t12.*t22.*t29.*t48.*(t2.*2.0+t17))./1.6e+1+(t12.*t29.*t37.*t47.*t51.*x3)./4.0;
end
end

function [d2,d3] = dI32_g(x2,x3,m2,m3,a2,a3,t,nu)
%dI32_g
%    [D2,D3] = dI32_g(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I32_g
t2 = a2+a3;
t3 = a2.*m2;
t4 = a2.*m3;
t5 = a3.*m2;
t6 = a3.*m3;
t7 = a2.*t;
t8 = a3.*t;
t9 = a2.*x2;
t10 = a2.*x3;
t11 = a3.*x2;
t12 = a3.*x3;
t13 = m3+x3;
t14 = a2.^2;
t15 = a2.^3;
t16 = a3.^2;
t17 = a3.^3;
t19 = m2.*2.0;
t20 = m3.*2.0;
t22 = m2.^2;
t23 = m3.^2;
t24 = t.^2;
t25 = x2.*2.0;
t26 = x3.*2.0;
t27 = x2.^2;
t28 = x3.^2;
t33 = 1.0./pi;
t34 = -a3;
t35 = 1.0./a2;
t36 = 1.0./a3;
t37 = nu-1.0;
t38 = -x2;
t40 = -x3;
t46 = m2.*x2.*-2.0;
t47 = m3.*x3.*-2.0;
t53 = m2./4.0;
t54 = m3./3.0;
t55 = m3./4.0;
t56 = x2./4.0;
t57 = x3.*(3.0./4.0);
t59 = m3.*x3.*(3.0./2.0);
t66 = nu-1.0./2.0;
t67 = nu-3.0./4.0;
t29 = t7.*2.0;
t30 = t8.*2.0;
t32 = t20.*x3;
t39 = -t25;
t41 = a2.*t13;
t42 = -t5;
t43 = -t9;
t44 = -t10;
t45 = -t12;
t48 = a2.*t4;
t49 = a2.*t7;
t50 = a3.*t8;
t51 = a3.*t12;
t52 = t13+t20;
t58 = a3.*t13.*2.0;
t60 = t8+t13;
t61 = a2+t34;
t62 = a2.*t36;
t63 = a3.*t35;
t64 = m2+t38;
t65 = m3+t40;
t68 = -t56;
t69 = t23./2.0;
t70 = 1.0./t37;
t73 = t54+x3;
t74 = t14+t16;
t76 = t7.^2;
t77 = t8.^2;
t85 = t4.*t13.*t14;
t99 = t55+t57;
t71 = t13.*t30;
t72 = -t48;
t75 = atan(t63);
t78 = t19+t39;
t79 = t64.^2;
t80 = a3.*t60;
t81 = a2.*t64;
t82 = a3.*t64;
t83 = a2.*t65;
t84 = m3.*t64;
t86 = t7.*t52;
t87 = t60.^2;
t88 = a3.*t65.*2.0;
t90 = t7+t64;
t91 = t8+t65;
t92 = 1.0./t60;
t94 = 1.0./t74.^2;
t95 = t34.*t64;
t98 = t29.*t64;
t100 = t24.*t74;
t101 = t13.*2.0+t30;
t103 = t8.*t17.*t64;
t105 = t7.*t73.*3.0;
t106 = t53+t68;
t108 = t7.*t64.*(3.0./2.0);
t110 = (t15.*t65)./4.0;
t113 = t4+t10+t11+t42;
t115 = a2.*t16.*t99;
t118 = t4+t11+t42+t44;
t138 = t3+t6+t12+t43+t49+t50;
t144 = t3+t6+t43+t45+t49+t50;
t89 = a3.*t81;
t93 = 1.0./t87;
t96 = -t83;
t97 = a2.*t78;
t102 = t90.^2;
t104 = a3.*t91;
t107 = t29+t78;
t109 = -t105;
t111 = (t14.*t82)./4.0;
t114 = t17.*t106;
t116 = t41+t95;
t119 = t113.^2;
t120 = t84+t86;
t121 = a2.*t113;
t122 = a3.*t113;
t123 = -t115;
t126 = t118.^2;
t136 = t49+t80+t81;
t142 = a2.*t138;
t143 = a3.*t138;
t146 = t138.^2;
t147 = t36.*t92.*t113;
t148 = t144.^2;
t150 = t28+t59+t69+t79+t108;
t168 = t23+t28+t32+t71+t76+t77+t79+t98;
t112 = -t111;
t117 = t58+t97;
t124 = t82+t96;
t125 = 1.0./t116;
t129 = t87+t102;
t130 = t88+t97;
t131 = t84+t109;
t132 = a3.*t14.*t120;
t133 = t51+t72+t89;
t140 = t49+t81+t104;
t141 = t92.*t116;
t149 = -t147;
t155 = a2.*t16.*t150.*2.0;
t160 = t121+t143;
t161 = t119+t146;
t167 = t126+t148;
t172 = 1.0./t168;
t128 = t.*t117;
t134 = t.*t130;
t135 = 1.0./t124;
t137 = 1.0./t129;
t139 = t17.*t131;
t145 = log(t141);
t151 = t62+t149;
t153 = t125.*t136;
t158 = -t155;
t165 = 1.0./t161;
t170 = 1.0./t167;
t173 = t172.^2;
t178 = t110+t112+t114+t123;
t152 = t151.^2;
t157 = atan(t153);
t164 = t135.*t140;
t169 = t22+t23+t27+t28+t32+t46+t100+t128;
t175 = t22+t23+t27+t28+t46+t47+t100+t134;
t154 = t152+1.0;
t162 = a2.*a3.*t157;
t166 = atan(t164);
t171 = log(t169);
t174 = 1.0./t169;
t176 = log(t175);
t177 = 1.0./t175;
t180 = t2.*t34.*t61.*t157.*t168;
t156 = log(t154);
t159 = 1.0./t154;
t181 = a2.*t16.*t168.*t171;
t182 = -t181;
t183 = t85+t103+t132+t139+t158+t180+t182;
d2 = t33.*t70.*(a2.*a3.*(a3.*t166-t124.*t170.*(a2.*t118+a3.*t144))+(a3.*t2.*t61.*t176)./4.0+(t2.*t61.*t107.*t124.*t177)./4.0).*(-1.0./4.0)-t33.*t36.*t66.*(t16+t16.*t145-(t16.*t156)./2.0+a2.*a3.*t75+a2.*t34.*t157-t74.*t87.*t137+a2.*t116.*t160.*t165+a3.*t141.*t151.*t159)-t33.*t67.*t70.*(t171.*(t17./4.0-(a3.*t14)./4.0)+a3.*(t162+t133.*t160.*t165)+t107.*t174.*t178)-(t33.*t70.*t94.*t172.*x3.*(-t6.*t16-t8.*t17+t34.*t48+a2.*t16.*(t7.*(3.0./2.0)+t78).*2.0+a2.*t16.*t107.*t171+a3.*t2.*t61.*t107.*t157+a2.*t16.*t107.*t168.*t174+a3.*t2.*t61.*t160.*t165.*t168))./4.0-(t33.*t70.*t94.*t107.*t173.*t183.*x3)./4.0;
if nargout > 1
    et1 = (t33.*t70.*(a2.*a3.*(a2.*t166-t124.*t170.*(a2.*t144+t34.*t118))+(a2.*t2.*t61.*t176)./4.0-(t2.*t61.*t124.*t177.*(t20-t26+t30))./4.0))./4.0+t33.*t67.*t70.*(-t171.*(t15./4.0+a2.*t16.*(3.0./4.0))+a3.*(t16.*t157+t133.*t165.*(t122-t142))+t101.*t174.*t178)-t33.*t36.*t66.*(t14.*t75-t14.*t157+t74.*atan(t90.*t92)+t80.*(a2.*t92-t93.*t116)+a2.*a3.*t145-(a2.*a3.*t156)./2.0-t60.*t74.*t90.*t137-a2.*t116.*t165.*(t122-t142)+a3.*t116.*t151.*t159.*(t62.*t92-t36.*t93.*t113))-(t33.*t70.*t94.*t172.*t183)./4.0;
    et2 = (t33.*t70.*t94.*t172.*x3.*(-t4.*t14+t7.*t17.*3.0+t7.*t14.*t34+a2.*t16.*(m3.*(3.0./2.0)+t26).*2.0+a2.*t16.*t101.*t171+a3.*t2.*t61.*t101.*t157+a2.*t16.*t101.*t168.*t174+a3.*t2.*t61.*t165.*t168.*(t122-t142)))./4.0+(t33.*t70.*t94.*t101.*t173.*t183.*x3)./4.0;
    d3 = et1+et2;
end
end

function [d2,d3] = dI32_h(x2,x3,m2,m3,a2,a3,t,nu)
%dI32_h
%    [D2,D3] = dI32_h(X2,X3,M2,M3,A2,A3,T,NU)


%x2- and x3-derivatives of I32_h
t2 = a2+a3;
t3 = a2.*m2;
t4 = a2.*m3;
t5 = a3.*m2;
t6 = a3.*m3;
t7 = a3.*t;
t8 = a2.*x2;
t9 = a2.*x3;
t10 = a3.*x2;
t11 = a3.*x3;
t12 = m3+x3;
t13 = a2.^2;
t14 = a2.^3;
t15 = a3.^2;
t16 = a3.^3;
t18 = m2.*2.0;
t19 = m3.*2.0;
t21 = m2.^2;
t22 = m3.^2;
t23 = t.^2;
t24 = x2.*2.0;
t25 = x3.*2.0;
t26 = x2.^2;
t27 = x3.^2;
t28 = a2.*t.*2.0;
t32 = 1.0./pi;
t33 = -a3;
t34 = nu-1.0;
t35 = -x2;
t37 = -x3;
t43 = m2.*x2.*-2.0;
t44 = m3.*x3.*-2.0;
t51 = m2./4.0;
t52 = m3./3.0;
t53 = m3./4.0;
t54 = x2./4.0;
t55 = x3.*(3.0./4.0);
t57 = m3.*x3.*(3.0./2.0);
t63 = nu-1.0./2.0;
t64 = nu-3.0./4.0;
t29 = t7.*2.0;
t31 = t19.*x3;
t36 = -t24;
t38 = a2.*t12;
t39 = -t5;
t40 = -t8;
t41 = -t9;
t42 = -t11;
t45 = a2.*t4;
t46 = t.*t13;
t47 = a3.*t7;
t48 = a3.*t11;
t49 = t12+t19;
t50 = t12.^2;
t56 = a3.*t12.*2.0;
t58 = t7+t12;
t59 = a2+t33;
t60 = m2+t35;
t61 = m3+t37;
t65 = -t54;
t67 = t22./2.0;
t68 = 1.0./t34;
t71 = t52+x3;
t72 = t13+t15;
t73 = t13.*t23;
t74 = t7.^2;
t82 = t4.*t12.*t13;
t93 = t53+t55;
t62 = 1.0./t50;
t66 = t+t60;
t69 = t12.*t29;
t70 = -t45;
t75 = t18+t36;
t76 = t60.^2;
t77 = a3.*t58;
t78 = a2.*t60;
t79 = a3.*t60;
t80 = a2.*t61;
t81 = m3.*t60;
t83 = a2.*t.*t49;
t84 = a3.*t61.*2.0;
t87 = t7+t61;
t88 = 1.0./t72.^2;
t89 = t33.*t60;
t92 = t28.*t60;
t94 = t23.*t72;
t95 = t12.*2.0+t29;
t96 = t7.*t16.*t60;
t98 = a2.*t.*t71.*3.0;
t99 = t51+t65;
t103 = (t14.*t61)./4.0;
t106 = t4+t9+t10+t39;
t109 = a2.*t15.*t93;
t112 = t4+t10+t39+t41;
t133 = t3+t6+t11+t40+t46+t47;
t138 = t3+t6+t40+t42+t46+t47;
t85 = t66.^2;
t86 = a3.*t78;
t90 = -t80;
t91 = a2.*t75;
t97 = a3.*t87;
t100 = t28+t75;
t101 = t.*t78.*(3.0./2.0);
t102 = -t98;
t104 = (t13.*t79)./4.0;
t108 = t16.*t99;
t110 = t38+t89;
t113 = t106.^2;
t114 = t81+t83;
t115 = a2.*t106;
t116 = a3.*t106;
t117 = -t109;
t121 = t112.^2;
t131 = t46+t77+t78;
t136 = a2.*t133;
t137 = a3.*t133;
t139 = t133.^2;
t140 = t138.^2;
t153 = t22+t27+t31+t69+t73+t74+t76+t92;
t105 = -t104;
t107 = t50+t85;
t111 = t56+t91;
t119 = t79+t90;
t120 = 1.0./t110;
t124 = t84+t91;
t125 = t22+t27+t31+t85;
t126 = t81+t102;
t127 = a3.*t13.*t114;
t128 = t48+t70+t86;
t135 = t46+t78+t97;
t141 = t27+t57+t67+t76+t101;
t146 = t115+t137;
t147 = t113+t139;
t152 = t121+t140;
t157 = 1.0./t153;
t118 = 1.0./t107;
t123 = t.*t111;
t129 = t.*t124;
t130 = 1.0./t119;
t132 = 1.0./t125;
t134 = t16.*t126;
t142 = t120.*t131;
t143 = a2.*t15.*t141.*2.0;
t150 = 1.0./t147;
t155 = 1.0./t152;
t158 = t157.^2;
t163 = t103+t105+t108+t117;
t144 = atan(t142);
t145 = -t143;
t149 = t130.*t135;
t154 = t21+t22+t26+t27+t31+t43+t94+t123;
t160 = t21+t22+t26+t27+t43+t44+t94+t129;
t151 = atan(t149);
t156 = log(t154);
t159 = 1.0./t154;
t161 = log(t160);
t162 = 1.0./t160;
t165 = t2.*t33.*t59.*t144.*t153;
t166 = a2.*t15.*t153.*t156;
t167 = -t166;
t168 = t82+t96+t127+t134+t145+t165+t167;
d2 = t32.*t68.*(a2.*a3.*(a3.*t151-t119.*t155.*(a2.*t112+a3.*t138))+(a3.*t2.*t59.*t161)./4.0+(t2.*t59.*t100.*t119.*t162)./4.0).*(-1.0./4.0)+t32.*t63.*(atan(t66./t12)+t12.*t66.*t118-(t12.*t132.*(t.*2.0+t75))./2.0)-t32.*t64.*t68.*(t156.*(t16./4.0-(a3.*t13)./4.0)+a3.*(a2.*a3.*t144+t128.*t146.*t150)+t100.*t159.*t163)-(t32.*t68.*t88.*t157.*x3.*(-t6.*t15-t7.*t16+t33.*t45+a2.*t15.*(t75+a2.*t.*(3.0./2.0)).*2.0+a2.*t15.*t100.*t156+a3.*t2.*t59.*t100.*t144+a2.*t15.*t100.*t153.*t159+a3.*t2.*t59.*t146.*t150.*t153))./4.0-(t32.*t68.*t88.*t100.*t158.*t168.*x3)./4.0;
if nargout > 1
    et1 = t32.*t63.*(log(t62.*t125)./2.0+t85.*t118-(t12.^3.*t132.*(1.0./t12.^3.*t125.*2.0-2.0./t12))./2.0)+(t32.*t68.*(a2.*a3.*(a2.*t151-t119.*t155.*(a2.*t138+t33.*t112))+(a2.*t2.*t59.*t161)./4.0-(t2.*t59.*t119.*t162.*(t19-t25+t29))./4.0))./4.0+t32.*t64.*t68.*(-t156.*(t14./4.0+a2.*t15.*(3.0./4.0))+a3.*(t15.*t144+t128.*t150.*(t116-t136))+t95.*t159.*t163)-(t32.*t68.*t88.*t157.*t168)./4.0;
    et2 = (t32.*t68.*t88.*t157.*x3.*(-t4.*t13-t7.*t14+a2.*t7.*t15.*3.0+a2.*t15.*(m3.*(3.0./2.0)+t25).*2.0+a2.*t15.*t95.*t156+a3.*t2.*t59.*t95.*t144+a2.*t15.*t95.*t153.*t159+a3.*t2.*t59.*t150.*t153.*(t116-t136)))./4.0+(t32.*t68.*t88.*t95.*t158.*t168.*x3)./4.0;
    d3 = et1+et2;
end
end

function [g1,g2,g3,g4,g5,g6,g7,g8] = dGreenParts(x2,x3,y2,y3,nu)
%dGreenParts
%    [G1,G2,G3,G4,G5,G6,G7,G8] = dGreenParts(X2,X3,Y2,Y3,NU)


%x2- then x3-derivatives of G22, G23, G32, G33 (point x, source y)
t2 = x3+y3;
t3 = nu.*2.0;
t4 = nu.*3.0;
t5 = nu.*4.0;
t7 = x2.*2.0;
t8 = x3.*2.0;
t9 = y2.*2.0;
t10 = y3.*2.0;
t13 = 1.0./pi;
t15 = nu-1.0;
t16 = -y2;
t18 = -y3;
t28 = nu-3.0./4.0;
t12 = t8.*y3;
t14 = -t4;
t17 = -t9;
t19 = -t10;
t20 = nu.*t3;
t21 = t10.*y3;
t22 = t2.^2;
t24 = t3-1.0;
t25 = t5-3.0;
t26 = t16+x2;
t27 = t18+x3;
t30 = t2.*y3.*-2.0;
t31 = 1.0./t15;
t32 = t7+t17;
t33 = t8+t19;
t34 = t26.^2;
t35 = t27.^2;
t36 = t22.*t25;
t37 = t14+t20+5.0./4.0;
t38 = t22+t34;
t39 = t25.*t34;
t40 = t34+t35;
t49 = t21+t30+t36;
t41 = 1.0./t38;
t44 = t12+t39;
t45 = 1.0./t40;
t42 = t41.^2;
t43 = t41.^3;
t46 = t45.^2;
t52 = (t26.*t45)./4.0;
t53 = (t27.*t45)./4.0;
t54 = t2.*t15.*t24.*t41;
t55 = (t25.*t26.*t41)./4.0;
t56 = (t25.*t27.*t41)./4.0;
t57 = t15.*t24.*t26.*t41;
t59 = (t28.*t32.*t45)./2.0;
t60 = (t28.*t33.*t45)./2.0;
t64 = t2.*t37.*t41;
t65 = (t32.*t37.*t41)./2.0;
t48 = t2.*t18.*t42.*x3;
t50 = t26.*t42.*x3.*y3;
t51 = t2.*t26.*t42.*y3;
t58 = t22.*t26.*t43.*x3.*y3.*4.0;
t61 = t2.*t12.*t26.*t32.*t43;
t63 = -t59;
g1 = (t13.*t31.*(t63+t65-(t32.*t35.*t46)./4.0+(t32.*t42.*t49)./4.0+t12.*t22.*t32.*t43))./2.0;
if nargout > 1
    t66 = (t26.*t27.*t32.*t46)./4.0;
    t67 = (t26.*t27.*t33.*t46)./4.0;
    t68 = (t2.*t25.*t26.*t27.*t42)./2.0;
    t69 = (t25.*t26.*t27.*t32.*t42)./4.0;
    g2 = t13.*t31.*(t48+t53+t54-t56+t61-t66+t69).*(-1.0./2.0);
end
if nargout > 2
    g3 = (t13.*t31.*(t48-t53+t54+t56+t61+t66-t69))./2.0;
end
if nargout > 3
    g4 = t13.*t31.*(t59-t65-(t32.*t45)./4.0+(t25.*t32.*t41)./4.0+(t32.*t34.*t46)./4.0-(t32.*t42.*t44)./4.0-t32.*t42.*x3.*y3+t32.*t34.*t43.*x3.*y3.*2.0).*(-1.0./2.0);
end
if nargout > 4
    g5 = (t13.*t31.*(-t60+t64+(t33.*t45)./4.0+(t41.*(t10-t2.*t25.*2.0))./4.0+t18.*t22.*t42+(t2.*t42.*t49)./2.0-(t33.*t35.*t46)./4.0+t30.*t42.*x3+t2.^3.*t43.*x3.*y3.*4.0))./2.0;
end
if nargout > 5
    t62 = -t58;
    g6 = (t13.*t31.*(t50+t51-t52+t55+t57+t62+t67-t68))./2.0;
end
if nargout > 6
    g7 = t13.*t31.*(t50+t51+t52-t55+t57+t62-t67+t68).*(-1.0./2.0);
end
if nargout > 7
    g8 = t13.*t31.*(t60-t64+(t41.*y3)./2.0-(t2.*t42.*t44)./2.0+t18.*t34.*t42+(t33.*t34.*t46)./4.0+t2.*t34.*t43.*x3.*y3.*4.0).*(-1.0./2.0);
end
end

function g=dGreen(x2,x3,y2,y3,nu)
% derivatives of the Green's function, as a cell array (see dGreenParts)
g=cell(1,8);
[g{:}]=dGreenParts(x2,x3,y2,y3,nu);
end
