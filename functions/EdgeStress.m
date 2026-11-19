function [Sxx,Sxz,Szz]=EdgeStress(m,x,z,nu,mu)
% EDGESTRESS  Stress due to an edge dislocation in a half-space.
%
%   [Sxx, Sxz, Szz] = EdgeStress(m, x, z, nu)
%   [Sxx, Sxz, Szz] = EdgeStress(m, x, z, nu, mu)
%
%   Computes the plane-strain stress at (x, z) due to a finite edge
%   dislocation (a dipping segment with uniform slip) in a homogeneous
%   elastic half-space, as the sum of two infinite edge dislocations of
%   opposite sign at its ends. All lengths must be in the same units;
%   otherwise the stresses are off by a constant.
%
% INPUTS
%   m       [1 x 5] dislocation: [depth of the updip end, horizontal
%           position of the updip end, length down dip, dip (degrees),
%           slip]
%   x, z    observation points (z is made negative inside the function)
%   nu      Poisson's ratio
%
% OPTIONS
%   mu      shear modulus; sets the units of the stress  [1]
%
% OUTPUTS
%   Sxx, Sxz, Szz  stress components, same size as x (units of mu)
%
% NOTES
%   - The vertical axis is negative down; the code takes -abs() of z and
%     of m(1), so their signs do not matter.
%
% SEE ALSO
%   computetractionkernels, LDstressHS
%
% AUTHORS
%   Not given in the source; added to I2EC by Rishav Mallick, 2023


if nargin <5	mu=1;end   
%Make sure vertical coordinates are negative	
	z=-abs(z);	Z=-abs(m(1));
   
%Define Constants	
	dip=m(4)*pi/180;	
	sd=sin(dip);	
	cd=cos(dip);	
	bx=m(5)*sd;	
	bz=m(5)*cd;	
	a1=[bx -bx]*mu/(2*pi*(1-nu));	
	a2=[bz -bz]*mu/(2*pi*(1-nu));	
	X=[m(2) m(2)+m(3)*cd];	
	Z=[Z Z-m(3)*sd];	
	Sxx=zeros(size(x));	
	Sxz=zeros(size(x));	
	Szz=zeros(size(x));   

%Compute stresses	
for i=1:2		
		zmZ=z-Z(i);
		zmZ2=zmZ.*zmZ;
		zpZ=z+Z(i);
		zpZ2=zpZ.*zpZ;
		xmX=x-X(i);
		xmX2=xmX.*xmX;		
		r12=zmZ2+xmX2;		
		r22=zpZ2+xmX2;		
		r122=r12.*r12;		
		r222=r22.*r22;		
		r223=r222.*r22;		
		b1=xmX.*(3*zmZ2+xmX2)./(r122);
		b2=xmX.*(3*zpZ2+xmX2)./(r222);	
		b3=4.0*Z(i)*z.*xmX.*(3.0*zpZ2-xmX2)./(r223);
		b4=zmZ.*(zmZ2-xmX2)./(r122);
		b5=zpZ.*(zpZ2-xmX2)./(r222);
		a3=(3.0*z+Z(i)).*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2-xmX2.*xmX2;
		b6=2.0*Z(i).*a3./(r223);

		c1=xmX.*(zmZ2-xmX2)./(r122);
		c2=xmX.*(zpZ2-xmX2)./(r222);
		a4=(2.0*Z(i)-z).*zpZ2+(3.0*z+2.0*Z(i)).*xmX2;
		c3=4.0*Z(i).*xmX.*a4./(r223);
		c4=zmZ.*(zmZ2+3.0*xmX2)./(r122);
		c5=zpZ.*(zpZ2+3.0*xmX2)./(r222);
		a5=zmZ.*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2+xmX2.*xmX2;
		c6=2.0*Z(i)*a5./(r223);	

		Sxx=Sxx+a1(i)*(c1-c2+c3)+a2(i)*(-c4+c5+c6);
		Sxz=Sxz+a1(i)*(b4-b5-c6)+a2(i)*(-c1+c2-b3);
		Szz=Szz+a1(i)*(-b1+b2+b3)+a2(i)*(-b4+b5-b6);	
end
