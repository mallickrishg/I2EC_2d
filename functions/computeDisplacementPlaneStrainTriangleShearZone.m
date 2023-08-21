function [u2,u3]=computeDisplacementPlaneStrainTriangleShearZone( ...
    x2,x3,A,B,C,e22,e23,e33,nu)
% function COMPUTEDISPLACEMENTPLANESTRAINTRIANGLESHEARZONE computes the
% displacement field associated with deforming triangle strain volume
% considering the following geometry using the tanh-sinh numerical
% quadrature.
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
% nu                 Poisson's ratio in the half space.
%
% Output:
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - March 7, 2018, Singapore.

assert(min(x3(:))>=0,'depth must be positive.');

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

% moment density
m22=lambda*ekk+2*e22;
m23=2*e23;
m33=lambda*ekk+2*e33;

u2= (m22*nC(1)+m23*nC(2))*T22(A,B) ...
   +(m23*nC(1)+m33*nC(2))*T32(A,B) ...
   +(m22*nA(1)+m23*nA(2))*T22(B,C) ...
   +(m23*nA(1)+m33*nA(2))*T32(B,C) ...
   +(m22*nB(1)+m23*nB(2))*T22(C,A) ...
   +(m23*nB(1)+m33*nB(2))*T32(C,A);

u3= (m22*nC(1)+m23*nC(2))*T23(A,B) ...
   +(m23*nC(1)+m33*nC(2))*T33(A,B) ...
   +(m22*nA(1)+m23*nA(2))*T23(B,C) ...
   +(m23*nA(1)+m33*nA(2))*T33(B,C) ...
   +(m22*nB(1)+m23*nB(2))*T23(C,A) ...
   +(m23*nB(1)+m33*nB(2))*T33(C,A);


pos=~isfinite(u2+u3);
if 0~=numel(x2(pos))
    [u2f,u3f]=computeDisplacementPlaneStrainTriangleShearZoneTanhSinh( ...
        x2(pos),x3(pos),A,B,C,e22,e23,e33,nu);
    u2(pos)=u2f;
    u3(pos)=u3f;

    fprintf('%d points corrected with double-exponential integration\n',numel(u2f));
end

       function u2=I22(A,B,t)
        
        % azimuthal vector
        a=(B-A)/2;
        a=a/norm(a);
        
        m2=(A(1)+B(1))/2;
        m3=(A(2)+B(2))/2;
        
        a2=a(1);
        a3=a(2);

        if 0==a2
            us=(-0.4e1.*(-x3+m3+t).*(nu-0.3e1./0.4e1).*log((m3.^2+(-2.*x3+2.*t).*m3+x3.^2-2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                +0.8e1.*((x2-m2).*atan(((-x3+m3+t)./(-x2+m2)))+t).*(nu-0.1e1./0.2e1))./pi./(-0.1e1+nu)./16;
            ui=(nu.^2-0.3e1./0.2e1.*nu+0.5e1./0.8e1).*((x3+m3+t).*log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                +(2.*m2-2.*x2).*atan(((x3+m3+t)./(-x2+m2)))-(2.*t))./pi./(-0.1e1+nu)./0.2e1;
            ui=ui+(x3.*(-x2+m2).*log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                +((-3+4.*nu).*m2.^2+(6-8.*nu).*x2.*m2+(-3+4.*nu).*x2.^2-2.*x3.^2).*atan(((x3+m3+t)./(-x2+m2))) ...
                -0.4e1.*(nu-0.3e1./0.4e1).*(-x2+m2).*t)./pi./(-1+nu)./(-x2+m2)./0.8e1;
            ui=ui-x3.*((-x2+m2).*(m2.^2-2.*x2.*m2+x2.^2+(x3+m3+t).^2).*log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                -x3.*(m2.^2-2.*x2.*m2+x2.^2+(x3+m3+t).^2).*atan(((x3+m3+t)./(-x2+m2))) ...
                +((-x2+m2).*(m2.^2-2.*x2.*m2+x2.^2+x3.*(x3+m3+t))))./pi./(-x2+m2)./(-1+nu)./(m2.^2-2.*x2.*m2+x2.^2+(x3+m3+t).^2)./0.4e1;
        else
        % source term
        us=-((t.*(nu-0.3e1./0.4e1).*a2.^4./0.2e1+(m2-x2).*(nu-0.3e1./0.4e1).*a2.^3./0.2e1 ...
            -(-0.2e1.*t.*(nu-0.3e1./0.4e1).*a3+(-0.5e1./0.4e1+nu).*(x3-m3)).*a3.*a2.^2./0.2e1 ...
            +(m2-x2).*(nu-0.1e1./0.4e1).*a3.^2.*a2./0.2e1-(nu-0.3e1./0.4e1).*a3.^3.*(x3-m3-a3.*t)./0.2e1) ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(-x3+m3)).*t+x2.^2 ...
                    -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2-0.2e1.*x3.*m3) ...
            +((nu-0.1e1./0.2e1).*a3.^2+a2.^2.*(-0.1e1+nu)).*((x3-m3).*a2+a3.*(m2-x2)) ...
            .*atan((-a2.^2.*t+(-m2+x2).*a2-(-x3+m3+a3.*t).*a3)./((-x3+m3).*a2-a3.*(m2-x2))) ...
            -t.*(a2.^2+a3.^2).*((nu-0.3e1./0.4e1).*a2.^2+(nu-0.1e1./0.2e1).*a3.^2)) ...
            ./(-0.1e1+nu)./pi./0.2e1;
        % first image term
        ui=((a2.^2.*t./0.2e1+(-x2./0.2e1+m2./0.2e1).*a2+a3.*(x3+m3+a3.*t)./0.2e1) ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
            -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2+0.2e1.*x3.*m3)+((x3+m3).*a2-a3.*(m2-x2)) ...
            .*atan((a2.^2.*t+(m2-x2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(m2-x2))) ...
            -(a2.^2+a3.^2).*t) ...
            .*(nu.^2-0.3e1./0.2e1.*nu+0.5e1./0.8e1)./(-0.1e1+nu)./pi;
        % second image term
        ui=ui+(-0.1e1./((-m2+x2).*a3+(x3+m3).*a2) ...
            .*(((-m2+x2).*a3+(x3+m3).*a2) ...
            .*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(m2-x2).*a2+x3.^2+2.*x3.*m3+m3.^2+(m2-x2).^2).*(a3.^3) ...
            .*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(x3+m3)).*t+x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2+2.*x3.*m3)) ...
            +(-a3.^4.*x3-3.*(m2-x2).*a2.*a3.^3+a2.^2.*(x3+3.*m3).*a3.^2-a2.^3.*(m2-x2).*a3+m3.*a2.^4) ...
            .*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(m2-x2).*a2+x3.^2+2.*x3.*m3+m3.^2+(m2-x2).^2) ...
            .*atan(((a2.^2.*t+(m2-x2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(m2-x2)))) ...
            -((-m2+x2).*a3+(x3+m3).*a2).*(-(t.*a3.^4.*x3)+((-3.*t.*(m2-x2).*a2-x3.^2-x3.*m3-(m2-x2).^2).*a3.^3) ...
            +0.3e1.*a2.*((t.*(x3+m3).*a2)-((m2-x2).*m3)./0.3e1) ...
            .*(a3.^2)+(a2.^2.*(t.*(m2-x2).*a2+x3.^2+3.*x3.*m3+2.*m3.^2+(m2-x2).^2).*a3)-(m3.*a2.^3.*(a2.*t+m2-x2)))).*x3 ...
            ./(-1+nu)./pi./(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(m2-x2).*a2+x3.^2+2.*x3.*m3+m3.^2+(m2-x2).^2)./0.4e1);
        % third image term
        ui=ui+(-0.4e1.*((-m2+x2).*a3+(x3+m3).*a2).*(-x3.*a3.^2./0.4e1 ...
            -a2.*(m2-x2).*(nu-0.3e1./0.4e1).*a3+a2.^2.*((x3+m3).*nu-x3-0.3e1./0.4e1.*m3)).*a3 ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
            -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2+0.2e1.*x3.*m3)+((0.4e1.*(m2-x2).^2.*nu-0.2e1.*x3.^2 ...
            -0.3e1.*(m2-x2).^2).*a3.^4-0.8e1.*a2.*((x3+m3).*nu-x3./0.2e1-0.3e1./0.4e1.*m3).*(m2-x2).*a3.^3 ...
            +0.4e1.*a2.^2.*((x3-m2+m3+x2).*(m2+x3+m3-x2).*nu ...
            -0.5e1./0.4e1.*x3.^2-x3.*m3+0.3e1./0.4e1.*(m2-m3-x2).*(m2+m3-x2)).*a3.^2 ...
            +0.8e1.*a2.^3.*(m2-x2).*((x3+m3).*nu-x3-0.3e1./0.4e1.*m3).*a3 ...
            -0.4e1.*a2.^4.*((x3+m3).^2.*nu-0.3e1./0.4e1.*x3.^2-0.2e1.*x3.*m3-0.3e1./0.4e1.*m3.^2)) ...
            .*atan((a2.^2.*t+(m2-x2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(m2-x2))) ...
            -0.4e1.*((-m2+x2).*a3+(x3+m3).*a2).*(nu-0.3e1./0.4e1).*a3.^2.*t) ...
            ./((-m2+x2).*a3+(x3+m3).*a2)./(-0.1e1+nu)./pi./0.8e1;
        end
        
        u2=us+ui;
    end

    function u2=T22(A,B)
        u2=(I22(A,B,norm(B-A)/2)-I22(A,B,-norm(B-A)/2));
    end

    function u3=I33(A,B,t)
        
        % azimuthal vector
        a=(B-A)/2;
        a=a/norm(a);
        
        m2=(A(1)+B(1))/2;
        m3=(A(2)+B(2))/2;
        
        a2=a(1);
        a3=a(2);
        
        if 0==a2
            us=(-0.4e1.*(-x3+m3+t).*(nu-0.3e1./0.4e1).*log((m3.^2+(-2.*x3+2.*t).*m3+x3.^2-2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                -0.8e1.*(-0.1e1+nu).*(-x2+m2).*atan(((-x3+m3+t)./(-x2+m2))) ...
                -(6.*t)+0.8e1.*t.*nu)./pi./(-0.1e1+nu)./0.16e2;
            
            % first image term
            ui=(nu.^2-0.3e1./0.2e1.*nu+0.5e1./0.8e1).*((x3+m3+t).*log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                +(2.*m2-2.*x2).*atan(((x3+m3+t)./(-x2+m2))) ...
                -(2.*t))./pi./(-0.1e1+nu)./0.2e1;
            
            % in the following expression the term x3.^2 disappears when
            % x2=m2 and the solution is finite.
            ui=ui+(-x3.*(-x2+m2).*log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)) ...
                -0.4e1.*((nu-0.3e1./0.4e1).*(m2.^2) ...
                        -0.2e1.*(nu-0.3e1./0.4e1).*x2.*m2 ...
                        +(nu-0.3e1./0.4e1).*(x2.^2)-(x3.^2)./0.2e1).*atan(((x3+m3+t)./(-x2+m2))) ...
                  )./pi./(-0.1e1+nu)./(-x2+m2)./0.8e1;
              
            % third image term
            ui=ui-x3.*((x3.^2+(2.*m3+2.*t).*x3+m2.^2-2.*x2.*m2+x2.^2+(m3+t).^2).*x3 ...
                        .*atan(((x3+m3+t)./(-x2+m2))) ...
                 +((x3.^2+(m3+t).*x3+(-x2+m2).^2).*(-x2+m2)))./pi ...
                 ./(x3.^2+(2.*m3+2.*t).*x3+m2.^2-2.*x2.*m2+x2.^2+(m3+t).^2)./(-1+nu)./(-x2+m2)./0.4e1;
        else
            
        % source term
        us=-((t.*(nu-0.3e1./0.4e1).*a2.^4./0.2e1 ...
                +(m2-x2).*(nu-0.3e1./0.4e1).*a2.^3./0.2e1 ...
                -a3.*(-0.2e1.*t.*(nu-0.3e1./0.4e1).*a3 ...
                +(nu-0.1e1./0.4e1).*(x3-m3)).*a2.^2./0.2e1 ...
                +(m2-x2).*(-0.5e1./0.4e1+nu).*a3.^2.*a2./0.2e1 ...
                -(nu-0.3e1./0.4e1).*a3.^3.*(x3-m3-a3.*t)./0.2e1) ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2 ...
                    +0.2e1.*a3.*(-x3+m3)).*t+x2.^2 ...
                    -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2 ...
                    -0.2e1.*x3.*m3) ...
            +((-0.1e1+nu).*a3.^2+a2.^2.*(nu-0.1e1./0.2e1)).*((x3-m3).*a2+a3.*(m2-x2)) ...
            .*atan((-a2.^2.*t+(-m2+x2).*a2-(-x3+m3+a3.*t).*a3)./((-x3+m3).*a2-a3.*(m2-x2))) ...
            -t.*(a2.^2.*(nu-0.1e1./0.2e1)+a3.^2.*(nu-0.3e1./0.4e1))) ...
            ./pi./(-0.1e1+nu)./0.2e1;
        
        % first image term
        ui=((nu.^2-0.3e1./0.2e1.*nu+0.5e1./0.8e1) ...
            .*((a2.^2.*t./0.2e1+(m2./0.2e1-x2./0.2e1).*a2+a3.*(x3+m3+a3.*t)./0.2e1) ...
            .*log(t.^2+((0.2e1.*m2-0.2e1.*x2).*a2 ...
                    +0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
                    -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2 ...
                    +0.2e1.*x3.*m3)+((x3+m3).*a2-a3.*(m2-x2)) ...
            .*atan((a2.^2.*t+(m2-x2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(m2-x2))) ...
            -(a2.^2+a3.^2).*t) ...
            ./pi./(-0.1e1+nu));
        % second image term
        ui=ui+((0.4e1.*a3.*((x3+m3).*a2-a3.*(-x2+m2)).*(((x3+m3).*nu-x3 ...
                    -0.3e1./0.4e1.*m3).*a2.^2 ...
                    -(-0.3e1./0.4e1+nu).*a3.*(-x2+m2).*a2-a3.^2.*x3./0.4e1) ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
                    -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2+0.2e1.*x3.*m3) ...
            +((0.4e1.*(x3+m3).^2.*nu-0.3e1.*m3.^2-0.8e1.*x3.*m3-0.3e1.*x3.^2).*a2.^4 ...
                -0.8e1.*a3.*(-x2+m2).*((x3+m3).*nu-x3-0.3e1./0.4e1.*m3).*a2.^3 ...
                -0.4e1.*((m2+x3+m3-x2).*(-m2+x3+m3+x2).*nu ...
                -0.5e1./0.4e1.*x3.^2-x3.*m3 ...
                +0.3e1./0.4e1.*(m2-m3-x2).*(m2+m3-x2)).*a3.^2.*a2.^2 ...
                +0.8e1.*a3.^3.*((x3+m3).*nu ...
                -0.3e1./0.4e1.*m3-x3./0.2e1).*(-x2+m2).*a2 ...
                +0.2e1.*(-0.2e1.*(-x2+m2).^2.*nu+x3.^2 ...
                +0.3e1./0.2e1.*(-x2+m2).^2).*a3.^4) ...
            .*atan((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2))) ...
            -0.4e1.*(-0.3e1./0.4e1+nu).*t.*((x3+m3).*a2-a3.*(-x2+m2)).*a2.^2) ...
            ./(-0.1e1+nu)./pi./((x3+m3).*a2-a3.*(-x2+m2))./0.8e1);
        % third image term
        ui=ui+(-0.3e1./0.4e1.*(-a3.*((-m2+x2).*a3+(x3+m3).*a2).*(a2.^2) ...
            .*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2 ...
                +2.*t.*(-x2+m2).*a2+x3.^2+2.*x3.*m3+m3.^2+(-x2+m2).^2) ...
            .*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(x3+m3)).*t ...
                    +x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2+2.*x3.*m3))./0.3e1 ...
            +(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+x3.^2+2.*x3.*m3+m3.^2+(-x2+m2).^2) ...
            .*((a3.^4.*x3)./0.3e1 ...
               -((-x2+m2).*a2.*a3.^3)./0.3e1 ...
               +(x3+m3./0.3e1).*(a2.^2).*(a3.^2)+(a2.^3.*(-x2+m2).*a3)./0.3e1 ...
               -(m3.*a2.^4)./0.3e1) ...
            .*atan(((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2)))) ...
            -(-(t.*a3.^4.*x3)+((-3.*t.*(-x2+m2).*a2-x3.^2-x3.*m3-(-x2+m2).^2).*a3.^3) ...
                +0.3e1.*((t.*(x3+m3).*a2)-((-x2+m2).*m3)./0.3e1).*a2.*(a3.^2) ...
                +((t.*(-x2+m2).*a2+x3.^2+3.*x3.*m3+2.*m3.^2+(-x2+m2).^2).*a2.^2.*a3) ...
                -(a2.^3.*m3.*(a2.*t+m2-x2))).*((-m2+x2).*a3+(x3+m3).*a2)./0.3e1) ...
            .*x3./((-m2+x2).*a3+(x3+m3).*a2) ...
            ./(-1+nu)./pi./(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+x3.^2+2.*x3.*m3+m3.^2+(-x2+m2).^2));
        end
        
        u3=us+ui;
    end

    function u3=T33(A,B)
        u3=(I33(A,B,norm(B-A)/2)-I33(A,B,-norm(B-A)/2));
    end



    function u3=I23(A,B,t)
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        
        m2=(A(1)+B(1))/2;
        m3=(A(2)+B(2))/2;
        
        a2=a(1);
        a3=a(2);
        
        if 0==a2
            us=-log((m3.^2+(-2.*x3+2.*t).*m3+x3.^2-2.*x3.*t+t.^2+(-x2+m2).^2)).*(-x2+m2)./pi./(-1+nu)./0.16e2;
            ui=-(nu-0.1e1./0.2e1).*((-m2./0.2e1+x2./0.2e1).*log(((m3.^2)+((2.*x3+2.*t).*m3)+(x3.^2)+(2.*x3.*t)+(t.^2)+(-x2+m2).^2)./((x3+m3+t).^2)) ...
                +(-x3-m3-t).*atan((-x2+m2)./(x3+m3+t)) ...
                +log((x2-m2)./(x3+m3+t)).*(-x2+m2))./pi;
            ui=ui+(log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)).*(-x2+m2) ...
                -0.4e1.*atan(((x3+m3+t)./(-x2+m2))).*x3).*(nu-0.3e1./0.4e1)./pi./(-0.1e1+nu)./0.4e1;
            ui=ui-x3.*((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2).*atan(((x3+m3+t)./(-x2+m2))) ...
                -((-x2+m2).*(m3+t)))./pi./(-1+nu) ...
                ./(m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)./0.4e1;
        else
        % source term
        us=(((x3-m3).*a2+a3.*(-x2+m2)).*(a2-a3).*(a2+a3) ...
            .*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(-x3+m3)).*t+x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2-2.*x3.*m3))./0.4e1 ...
            +a2.*(((x3-m3).*a2+a3.*(-x2+m2)) ...
            .*atan(((-t.*a2.^2+(-m2+x2).*a2-a3.*(a3.*t-x3+m3))./((-x3+m3).*a2-a3.*(-x2+m2)))) ...
            -t./0.2e1).*a3)./(-1+nu)./pi./0.4e1;
        
        if 0~=a3 
        % first image term
        ui=(-(a3.*(-m2+x2)+(x3+m3).*a2).*a3 ...
            .*log(((-a2./a3+(a2.*x3+a2.*m3+a3.*x2-a3.*m2)./a3./(x3+m3+a3.*t)).^2+1))./0.2e1 ...
            -(a3.*(-m2+x2)+(x3+m3).*a2).*a2.*atan(((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2))))+(a3.*(-m2+x2)+(x3+m3).*a2).*a3 ...
            .*log(((a3.*(-m2+x2)+(x3+m3).*a2)./(x3+m3+a3.*t)))+(a3.*(-m2+x2)+(x3+m3).*a2).*a2 ...
            .*atan((a3./a2))+(a2.^2+a3.^2).*(x3+m3+a3.*t) ...
            .*atan(((-x2+m2+a2.*t)./(x3+m3+a3.*t)))).*(-0.1e1./0.2e1+nu)./pi./(a2.^2+a3.^2)./a3;
        else
           ui=-1/2.*((x3+m3).*log(((-x2+m2+t).^2+x3.^2+2.*x3.*m3+m3.^2)./(x3+m3).^2)-2*atan((-x2+m2+t)./(x3+m3)).*(-x2+m2+t)).*(nu-1/2)/pi;
        end
        % second image term
        ui=ui+((-0.3e1./0.4e1+nu).*(((-x2./0.4e1+m2./0.4e1).*a3.^3-0.3e1./0.4e1.*(x3+m3./0.3e1).*a2.*a3.^2-a2.^2.*(-x2+m2).*a3./0.4e1-a2.^3.*(x3-m3)./0.4e1) ...
            .*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
            -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2+0.2e1.*x3.*m3)+((a3.^2.*x3+(-x2+m2).*a2.*a3-a2.^2.*m3) ...
            .*atan((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2))) ...
            +a2.*t.*(a2.^2+a3.^2)./0.2e1).*a3)./(-0.1e1+nu)./pi);
        % third image term
        ui=ui+(-0.1e1./((a2.^2+a3.^2).^2).*x3.*(a2.*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2).*(a3.^2) ...
            .*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(x3+m3)).*t+x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2+2.*x3.*m3)) ...
            +(a2-a3).*(a2+a3).*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2).*a3 ...
            .*atan(((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2)))) ...
            -(t.*(-x2+m2).*a3.^4)+(0.3e1.*t.*(x3+m3./0.3e1).*a2-(m3.*(-x2+m2))).*(a3.^3) ...
            +0.2e1.*a2.*(0.3e1./0.2e1.*t.*(-x2+m2).*a2+(m3.^2)./0.2e1 ...
            +0.3e1./0.2e1.*x3.*m3+(x3.^2)+((-x2+m2).^2)).*(a3.^2)-(a2.^2.*(t.*(3.*m3+x3).*a2+m3.*(-x2+m2)).*a3)-(m3.*a2.^3.*(x3+m3))) ...
            ./(-1+nu)./pi./(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2)./0.4e1);
        end
        
        u3=us+ui;
    end

    function u3=T23(A,B)
        u3=(I23(A,B,norm(B-A)/2)-I23(A,B,-norm(B-A)/2));
    end

    function u2=I32(A,B,t)
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        
        m2=(A(1)+B(1))/2;
        m3=(A(2)+B(2))/2;
        
        a2=a(1);
        a3=a(2);
        
        if 0==a2
            us=-log((m3.^2+(-2.*x3+2.*t).*m3+x3.^2-2.*x3.*t+t.^2+(-x2+m2).^2)).*(-x2+m2)./pi./(-1+nu)./0.16e2;
            ui=(nu-0.1e1./0.2e1).*((-m2./0.2e1+x2./0.2e1).*log(((m3.^2)+((2.*x3+2.*t).*m3)+(x3.^2)+(2.*x3.*t)+(t.^2)+(-x2+m2).^2)./((x3+m3+t).^2)) ...
                +(-x3-m3-t).*atan((-x2+m2)./(x3+m3+t)) ...
                +log((x2-m2)./(x3+m3+t)).*(-x2+m2))./pi;
            ui=ui+(log((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)).*(-x2+m2) ...
                -0.4e1.*atan(((x3+m3+t)./(-x2+m2))).*x3).*(nu-0.3e1./0.4e1)./pi./(-0.1e1+nu)./0.4e1;
            ui=ui+x3.*((m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2).*atan(((x3+m3+t)./(-x2+m2)))-((-x2+m2).*(m3+t)))./pi./ ...
                (m3.^2+(2.*x3+2.*t).*m3+x3.^2+2.*x3.*t+t.^2+(-x2+m2).^2)./(-1+nu)./0.4e1;
        else
        %sourceterm
        us=(((x3-m3).*a2+a3.*(-x2+m2)).*(a2-a3).*(a2+a3).*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(-x3+m3)).*t+x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2-2.*x3.*m3))./0.4e1 ...
            +a2.*(((x3-m3).*a2+a3.*(-x2+m2)).*atan(((-t.*a2.^2+(-m2+x2).*a2-a3.*(a3.*t-x3+m3))./((-x3+m3).*a2-a3.*(-x2+m2))))-((a2.^2+a3.^2).*t)./0.2e1).*a3)./(-1+nu)./pi./0.4e1;
        
        % first image term
        if 0~=a3
        ui=-(-(a3.*(-m2+x2)+(x3+m3).*a2).*a3.*log(((-a2./a3+(a2.*x3+a2.*m3+a3.*x2-a3.*m2)./a3./(x3+m3+a3.*t)).^2+1))./0.2e1 ...
             -(a3.*(-m2+x2)+(x3+m3).*a2).*a2.*atan(((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2)))) ...
             +(a3.*(-m2+x2)+(x3+m3).*a2).*a3.*log(((a3.*(-m2+x2)+(x3+m3).*a2)./(x3+m3+a3.*t)))+(a3.*(-m2+x2)+(x3+m3).*a2).*a2 ...
             .*atan((a3./a2))+(a2.^2+a3.^2).*(x3+m3+a3.*t).*atan(((-x2+m2+a2.*t)./(x3+m3+a3.*t)))).*(-0.1e1./0.2e1+nu)./pi./a3;
        else
                ui=1/2*((x3+m3).*log(((-x2+m2+t).^2+x3.^2+2*x3.*m3+m3.^2)./(x3+m3).^2)-2.*atan((-x2+m2+t)./(x3+m3)).*(-x2+m2+t)).*(nu-1/2)/pi;
        end
        
        %second image term
        ui=ui+((-0.3e1./0.4e1+nu).*(((-x2./0.4e1+m2./0.4e1).*a3.^3-0.3e1./0.4e1.*(x3+m3./0.3e1).*a2.*a3.^2-a2.^2.*(-x2+m2).*a3./0.4e1 ...
              -a2.^3.*(x3-m3)./0.4e1).*log((a2.^2+a3.^2).*t.^2+((0.2e1.*m2-0.2e1.*x2).*a2+0.2e1.*a3.*(x3+m3)).*t+x2.^2 ...
              -0.2e1.*x2.*m2+x3.^2+m2.^2+m3.^2+0.2e1.*x3.*m3)+((a3.^2.*x3+(-x2+m2).*a2.*a3-a2.^2.*m3) ...
              .*atan((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2)))+a2.*t.*(a2.^2+a3.^2)./0.2e1).*a3)./(-0.1e1+nu)./pi);
          
        %third image term
        ui=ui-(-0.1e1./((a2.^2+a3.^2).^2).*x3.*(a2.*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2).*(a3.^2) ...
            .*log(((a2.^2+a3.^2).*t.^2+((2.*m2-2.*x2).*a2+2.*a3.*(x3+m3)).*t+x2.^2-2.*x2.*m2+x3.^2+m2.^2+m3.^2+2.*x3.*m3)) ...
            +(a2-a3).*(a2+a3).*(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2).*a3 ...
            .*atan(((t.*a2.^2+(-x2+m2).*a2+a3.*(x3+m3+a3.*t))./((x3+m3).*a2-a3.*(-x2+m2)))) ...
            -(t.*(-x2+m2).*a3.^4) ...
            +(0.3e1.*t.*(x3+m3./0.3e1).*a2-(m3.*(-x2+m2))).*(a3.^3) ...
            +0.2e1.*a2.*(0.3e1./0.2e1.*t.*(-x2+m2).*a2+(m3.^2)./0.2e1 ...
            +0.3e1./0.2e1.*x3.*m3+(x3.^2)+((-x2+m2).^2)).*(a3.^2)-(a2.^2.*(t.*(3.*m3+x3).*a2+m3.*(-x2+m2)).*a3)-(m3.*a2.^3.*(x3+m3))) ...
            ./(-1+nu)./pi ...
            ./(a3.^2.*t.^2+2.*t.*(x3+m3).*a3+a2.^2.*t.^2+2.*t.*(-x2+m2).*a2+m3.^2+2.*x3.*m3+x3.^2+(-x2+m2).^2)./0.4e1);
        end
        u2=us+ui;
    end

    function u2=T32(A,B)
        u2=(I32(A,B,norm(B-A)/2)-I32(A,B,-norm(B-A)/2));
    end

end
