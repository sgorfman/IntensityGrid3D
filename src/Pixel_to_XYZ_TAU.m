function [HX,HY,HZ,TTH,Psi,TTH_hor,Zrel,Oblique,Rx,Ry,Rz] = Pixel_to_XYZ_TAU(X,Y,P)
%This function converts the X,Y coordinates of the PILATUS to the
%reciprocal space coorditanes HX, HY, HZ 
% 
% The coordinate system is such that X-axis looks to the X-ray
%source, Y --> X is the same as the positive TTH rotation of the
%detector

%P are parameter vectors of the PILATUS, which includes
%P(1) - distance to the detector in mm 
%P(2:3) - Position of the beam in pixel when detector angle (TTH) is zero
%P(4) - Pixel size in mm
%P(5) - detector angle, current value of TTH in degrees
%P(6) - wavelength in Angstrom

Roll = 0; 
if (numel(P) == 10)
    Roll = P(10);
end

d0 = P(1); %The distance from the rotation center to the beam
alpha = P(5); %The angle of the detector

D0 = [-d0*cosd(alpha),d0*sind(alpha),0]; %Cartesian coordinates of the detector center (The point on the detector where the beam would hit)
D0u = [-cosd(alpha),sind(alpha),0]; %The unit vector pointing at the detector center

Yd = [0,0,1];  %The orientation of the detector Yd axis
Xd = [sind(alpha),cosd(alpha),0]; %The orientation of the detector Xd axis
%Now let us correct this for the detector roll around the detector
M = [cosd(Roll),-sind(Roll); sind(Roll),cosd(Roll)];
MXY = [Xd',Yd']*M;
Xd = MXY(:,1)';
Yd = MXY(:,2)';
%================================================================================

Rx = D0(1) + P(4)*(X - P(2))*Xd(1) + P(4)*(Y - P(3))*Yd(1);
Ry = D0(2) + P(4)*(X - P(2))*Xd(2) + P(4)*(Y - P(3))*Yd(2);
Rz = D0(3) + P(4)*(X - P(2))*Xd(3) + P(4)*(Y - P(3))*Yd(3);

%Rx, Ry, Rz are the coordinates of the Bragg peak on the detector in the
%XYZ coordinate system

%Let us add one line which takes into account the fact that the crystal may
%be displaced from the center of rotation and the displacement is
%P(7),P(8),P(9)

%This option is going to be implemented only if the number of columns of
%the vector P is greater than 6
if (size(P,2)>6)
    Rx = Rx - P(7);
    Ry = Ry - P(8);
    Rz = Rz - P(9);
end
%Otherwise this correction is not going to be taken into account. 

L = (Rx.^2 + Ry.^2 + Rz.^2).^(-0.5); 

Rx = Rx.*L; 
Ry = Ry.*L;
Rz = Rz.*L;

Oblique = acosd(Rx*D0u(1) + Ry*D0u(2) + Rz*D0u(3)); %This is the angle between the diffracted X-ray beam and the normal to the detector

TTH = acosd(-Rx); %This is the scattering angle 
TTH_hor = acosd(-Rx.*((Rx.^2 + Ry.^2).^(-0.5))); 
Zrel = (Y - P(3)); 

HX = (Rx + 1) / P(6);
HY =  Ry / P(6);
HZ =  Rz / P(6);
Psi = atan2d(HZ,HY); %The asimuthal angle for the integration

%HX, HY and HZ are the orientation of the scattering vector in the
%Cartesian coordinate system
