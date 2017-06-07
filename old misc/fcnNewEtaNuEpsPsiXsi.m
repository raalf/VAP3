function [normal, nu, epsilon, psi, xsi, phiLE, phiMID, phiTE, eta] = fcnNewEtaNuEpsPsiXsi(Temp, US_DVE, xo, normal, phiLE, phiMID, phiTE)
% //This function computes the roll, pitch and sweep angles, as well the
% //halfspan of a wake DVE after it has been "relaxed".
% //nu and epsilon are computed from the normal of the DVE plane.  The normal
% //is the crossproduct of the vector delX1 and delX2.  delX1 is the vector
% //from the half point of the trailing edge of the next upstream DVE to the new
% //reference point. delX2 is the vector between the left side-edge point and the
% //reference point.  This vector also is used to determine the change in sweep
% //and the new halfspan of the DVE.
% //Yaw, psi, is kept constant as is the halfchord xsi.
% //
% //
% //input:
% //	wakeDVE	- DVE of interest
% //	US_DVE	- DVE upstream
% //
% //output:
% //updated values for wakeDVE
% //	nu		- roll angle
% //	epsilon - pitch (incidence) angle
% //	psi		- yaw angle
% //	xsi		- half chord length (in local, psi-rotated system)
% //	eta		- half span


%% normal
tempA(1) = US_DVE.xsi;
% tempA(1) = US_DVE.xsi;
tempA(2) = 0;
tempA(3) = 0;
% Trailing edge midpoint of upstream DVE

xteUS = star_glob(tempA, US_DVE.nu, US_DVE.epsilon, US_DVE.psi);

xteUS = xteUS + US_DVE.xo;

delX1 = xo - xteUS;
delX2 = normal;

tempA = cross(delX1, delX2);

normal = tempA./norm(tempA);

%% nu
if abs(normal(3)) > Temp.DBL_EPS
    nu = -atan2d(normal(2),normal(3));
%     if normal(3) < 0
% %         nu = nu + pi;
%     nu = nu + 180;
%     end
% else
%     if normal(2) > 0
% %         nu = 0.5*pi;
%         nu = -90;
%     else
% %         nu = 0.5*pi;
%         nu = 90;
%     end
end

%% epsilon
epsilon = asind(normal(1));
% epsilon = asin(normal(1));

%% psi

% delXSI1 = glob_star(delX1, US_DVE.nu, US_DVE.epsilon, 0);
delXSI1 = glob_star(delX1, nu, epsilon, 0);

if delXSI1(1)*delXSI1(1) > Temp.DBL_EPS
%     psi = atand(delXSI1(2)/delXSI1(1));
    psi = atan2d(delXSI1(2),delXSI1(1));
%     if delXSI1(1) < 0
% %         psi = psi + pi;
%         psi = psi + 180;
%     end
else
%     psi = 0.5*pi*abs(delXSI1(2))/delXSI1(1);
    psi = rad2deg(0.5*pi*abs(delXSI1(2))/delXSI1(1));
end

%% xsi
xsi = sqrt(delXSI1(1)*delXSI1(1) + delXSI1(2)*delXSI1(2));

%% sweep

delXSI2 = glob_star(delX2, nu, epsilon, psi);

if delXSI2(2) < 0
   disp('Problem in line 74 of fcnNewEtaNuEpsPsiXsi (relaxed wake stuff)')
   % See line 805-ish of wake_geometry.cpp by Bramesfeld
end

delPhi = atand(delXSI2(1)/delXSI2(2)) - phiMID;

phiLE = phiLE + delPhi;
phiMID = phiMID + delPhi;
phiTE = phiTE + delPhi;

%% eta
eta = delXSI2(2);


end

