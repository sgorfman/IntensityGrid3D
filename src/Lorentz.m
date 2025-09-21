function [Vp,L,M1,M2] = Lorentz(X,Y,P,Omega,CHI,Phi,dPhi,dOmega)
dS = P(4)^2; %The area of one pixel in mm^2
dAngle = max([dPhi,dOmega]); %In fact we assume that one of these angles is equal to zero, becasue only one angle is changed during the diffractometer rotation
M1 = RotationMatrix_TAU([Omega,CHI,Phi]); 
M2 = RotationMatrix_TAU([Omega+dOmega,CHI,Phi+dPhi]); 
M1 = inv(M1);
M2 = inv(M2);
VXYphi = dS*dAngle;

[HX,HY,HZ] = Pixel_to_XYZ_TAU(X,Y,P); 
HX1 = M1(1,1)*HX + M1(1,2)*HY + M1(1,3)*HZ;
HY1 = M1(2,1)*HX + M1(2,2)*HY + M1(2,3)*HZ;
HZ1 = M1(3,1)*HX + M1(3,2)*HY + M1(3,3)*HZ;

HX2 = M2(1,1)*HX + M2(1,2)*HY + M2(1,3)*HZ;
HY2 = M2(2,1)*HX + M2(2,2)*HY + M2(2,3)*HZ;
HZ2 = M2(3,1)*HX + M2(3,2)*HY + M2(3,3)*HZ;
    
[DHX_X,DHX_Y] = gradient(HX1,1,1);  DHX_Phi = (HX2 - HX1);
[DHY_X,DHY_Y] = gradient(HY1,1,1);  DHY_Phi = (HY2 - HY1);
[DHZ_X,DHZ_Y] = gradient(HZ1,1,1);  DHZ_Phi = (HZ2 - HZ1);

Vp = DHX_X.*(DHY_Y.*DHZ_Phi - DHY_Phi.*DHZ_Y) - DHX_Y.*(DHY_X.*DHZ_Phi - DHY_Phi.*DHZ_X) + DHX_Phi.*(DHY_X.*DHZ_Y - DHY_Y.*DHZ_X);
Vp = abs(Vp); %This is reciprocal space volume itself
L = VXYphi*Vp.^(-1);  %This is the Lorentz-factor that describes relationship between dX dY dPhi and the reciprocal space volume

M1 = inv(M1);
M2 = inv(M2);