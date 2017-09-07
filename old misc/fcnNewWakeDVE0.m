function [normal, nu, epsilon, psi, phiLE, phiMID, phiTE, eta, xo, u] = fcnNewWakeDVE0(Temp, US_DVE, xsi)

% 	//This function computes the reference point, angles, span and left
% 	//edge points of the first row of wake DVEs (time index = 0).
% 	//These are aligned with the free stream and are attached at their
% 	//leading edge to the trailing edges of the the upstream DVEs (index 1)
% 	//added June, 11 2005, G.B.
% 	//
% 	//input:
% 	//	info		- general information
% 	//	wakeDVE0	- first row of wake DVEs, time index = 0
% 	//	wakeDVE1	- row of wake DVEs directly upstream of wakeDVE0
% 	//
% 	//ouptut:
% 	//  wakeDVE0	- updated xo, nu, epsilon, psi, eta, xleft

tempA(1) = US_DVE.xsi;
% tempA(1) = US_DVE.xsi;
tempA(2) = 0;
tempA(3) = 0;
% Trailing edge midpoint of upstream DVE

xteUS = star_glob(tempA, US_DVE.nu, US_DVE.epsilon, US_DVE.psi);

xteUS = xteUS + US_DVE.xo;
% //the new ref. point is located along the extension of the free stream
xo = xteUS + (xsi.*Temp.u);

% 		//determining sweep angles and half span
%
% 		//computing vector along trailing edge of wing, xte
% 		//in local reference frame
tempA(1) = tand(US_DVE.phiTE)*US_DVE.eta;
tempA(2) = US_DVE.eta;
tempA(3) = 0;
% //transformation into global ref. frame
xte = star_glob(tempA, US_DVE.nu, US_DVE.epsilon, US_DVE.psi);

ABS_xte = norm(xte);

% 		//computing leading edge sweep, phiLE
% 		//phiLE is the angle between the leading edge, vector xte,
% 		//and the xsi_star axis, or the vector u

tempS = dot(xte, Temp.u)/ABS_xte;

phiLE = asind(tempS);
% 		//the remaining angles are of equal value since
% 		//they don't really matter for semi-infinite sheet
phiMID = phiLE;
phiTE = phiLE;

% //DVE half-span, projection of xte onto eta_star
eta = ABS_xte*cosd(phiLE);

% 		//determining nu, epsilon, and psi
% 
% 		//computing the normal of DVE,

tempA = cross(Temp.u, xte);
normal = tempA./norm(tempA);

if abs(normal(3)) > Temp.DBL_EPS
    nu = -atand(normal(2)/normal(3));
% nu = -atan(normal(2)/normal(3));
    if normal(3) < 0
%         nu = nu + pi;
    nu = nu + 180;
    end
else
    if normal(2) > 0
%         nu = 0.5*pi;
        nu = -90;
    else
%         nu = 0.5*pi;
        nu = 90;
    end
end

epsilon = asind(normal(1));

tempS = (Temp.u(2)*cosd(nu) + Temp.u(3)*sind(nu));
psi = asind(tempS);

u = US_DVE.u;




end

