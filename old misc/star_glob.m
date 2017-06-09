function [ x ] = star_glob( xsi,nu,eps,psi )
%STAR_GLOB
% nu = deg2rad(nu);
% eps = deg2rad(eps);
% psi = deg2rad(psi);
% double cnu = cos(nu), snu=sin(nu);
% double ceps = cos(eps), seps=sin(eps);
% double cpsi = cos(psi), spsi=sin(psi);
% double rot[3][3];		//rotation matrix

% //the rotation matrix:
% 	rot[0][0] = cpsi*ceps;
% 				rot[0][1] = cpsi*seps*snu+spsi*cnu;
% 							rot[0][2] =-cpsi*seps*cnu+spsi*snu;
% 	rot[1][0] =-spsi*ceps;
% 				rot[1][1] =-spsi*seps*snu+cpsi*cnu;
%  							rot[1][2] = spsi*seps*cnu+cpsi*snu;
% 	rot[2][0] = 	 seps;
% 				rot[2][1] = 	-ceps*snu;
% 			 					rot[2][2] =		 ceps*cnu;
% 
% 	//transforming x into xsi*
% 	x[0] = xsi[0]*rot[0][0] + xsi[1]*rot[1][0] + xsi[2]*rot[2][0];
% 	x[1] = xsi[0]*rot[0][1] + xsi[1]*rot[1][1] + xsi[2]*rot[2][1];
% 	x[2] = xsi[0]*rot[0][2] + xsi[1]*rot[1][2] + xsi[2]*rot[2][2];

cnu = cosd(nu);
snu= sind(nu);
ceps = cosd(eps);
seps= sind(eps);
cpsi = cosd(psi);
spsi= sind(psi);

rot(1,1) = cpsi*ceps;
rot(1,2) = cpsi*seps*snu+spsi*cnu;
rot(1,3) =-cpsi*seps*cnu+spsi*snu;
rot(2,1) =-spsi*ceps;
rot(2,2) =-spsi*seps*snu+cpsi*cnu;
rot(2,3) = spsi*seps*cnu+cpsi*snu;
rot(3,1) = 	 seps;
rot(3,2) = 	-ceps*snu;
rot(3,3) =	ceps*cnu;

x(1) = xsi(1)*rot(1,1) + xsi(2)*rot(2,1) + xsi(3)*rot(3,1);
x(2) = xsi(1)*rot(1,2) + xsi(2)*rot(2,2) + xsi(3)*rot(3,2);
x(3) = xsi(1)*rot(1,3) + xsi(2)*rot(2,3) + xsi(3)*rot(3,3);
end

