function [a3x, b3x, c3x] = fcnDVEInduction(Temp, xA, xo, nu, eps, phiLE, phiTE, ...
    psi, eta, xsi, DVE_type, singfct)
% Facsimile of DVEInduction from induced_velocity.cpp

%% Comments from FreeWake:
% //this function computes influence coefficients in point i due to
% //a distributed vorticity element.
% //The element has two vortex lines that are apart by 2*xsi and have a span
% //of 2*eta.  A vortex sheet is inbetween the two vortex lines.  If the
% //leading vortex line has a circulation of A+B*eta+C*eta^2 the aft line's
% //circulation is the negative value.  The sheet's vorticity is the derivative
% //B+2*eta*C.
% //The influence coefficient of a DVE is the combination of four separate
% //influence coefficients computed according to KHH: two of lifting lines
% //and two of semi-infinite vortex sheets, one starting at the leading edge
% //of the DVE and one of negative magnitude starting at the trailing edge of
% //the DVE.
% //the return values are in global coordinates
% //
% //Input:
% //	xA		  	point where velocity is induced in global coordinates
% //  xo 	  		reference point of DVE that induces
% //  nu,eps		DVE roll, pitch angles
% //  phiLE,phiTE	DVE leading and trailing edge sweeps
% //	psi			DVE yaw angle
% //  eta,xsi		DVE half-span and half-chord
% //  DVE_type 	type of DVE,
% //			DVE_type = 0 DVE has vortex filaments at leading and
% //						 trailing edge, usually lifting surface DVE
% //			DVE_type = 1 DVE has no vortex filaments at leading and
% //						 trailing edge, usually wake DVE
% //			DVE_type = 2 DVE has vortex filament at leading edge,
% //						 but not at trailing edge
% //			DVE_type =-2 DVE has vortex filament at trailing edge,
% //						 but not at leading edge
% //			DVE_type = 3 DVE is a semi infinite vortex sheet without a
% //						 vortex filaments at leading and trailing edge
% //			DVE_type =-3 DVE is a semi infinite vortex sheet with a
% //						 vortex filaments at its leading edge
% //			DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
% //						 xsi aft of the ref. pt. (for CD computation along TE)
% //			DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
% //						 also vortex filament at LE (CD computation along TE)
% //  singfct rate at which singularity at edge of vortex sheet decays
% //			## singfct added 2/8/05 G.B.  This factor is only needed when
% //			relaxing the wake, otherwise it is set to zero (0)
% //
% //Output:
% //	a3x[3]	induced coefficient of DVE in global system
% //	b3x[3]	induced coefficient of DVE in global system
% //	c3x[3]	induced coefficient of DVE in global system

%% Function Body

% 	//Expressing point A with respect to the DVE reference point in the
% 	//local DVE coordinates
% 	//1. xA-Xo

rA = xA - xo;

% fprintf('rA: %f %f %f\n', rA(1), rA(2), rA(3));

% //rotation to loca DVE reference frame, first nu, then epsilon, then psi
xsiA = glob_star(rA,nu,eps,psi);

% fprintf('rA: %f %f %f xsiA: %f %f %f\n', rA(1), rA(2), rA(3), xsiA(1), xsiA(2), xsiA(3));
% fprintf('rA:\t%f\t%f\t%f\tnu:\t%f\teps:\t%f\tpsi:\t%f\txsiA:\t%f\t%f\t%f\n',rA(1), rA(2), rA(3),nu,eps,psi, xsiA(1), xsiA(2), xsiA(3));
% 	//computing the leading edge influence of the DVE //

% 	//vortex system reference point midspan of DVE leading edge,
% 	//unless when needed for CDiEppler computation. added 8/13/05 G.B.
if DVE_type == 4
    xsio(1) = 0.5*xsi;
else
    xsio(1) = -xsi;
end

xsio(2) = 0; xsio(3) = 0;

% 	//computes influence coefficients, a1,b1, and c1, of
% 	//bound leading edge vortex of element on point xsiA
% 	//nu=0 since rotation already done

a1le(1) = 0;  a1le(2) = 0; a1le(3) = 0;
b1le(1) = 0;  b1le(2) = 0; b1le(3) = 0;
c1le(1) = 0;  c1le(2) = 0; c1le(3) = 0;
    
if DVE_type == 0 || DVE_type == 2 || DVE_type == -3 || DVE_type == -4
    [a1le, b1le, c1le] = fcnBoundVortexInduction(Temp, xsiA, xsio, 0, phiLE, eta);
end

% 	//computes influence coefficients, a2,b2, and c2, of
% 	//semi-infinite vortex sheet starting at leading edge of element
% 	//on point xsiA
% 	//nu=0 since rotation already done
%phile goes in as deg.
[a2le, b2le, c2le] = fcnVortexSheetInduction(Temp, xsiA, xsio, 0, phiLE, eta, singfct);

a3xi = a1le;
% disp('b1le: ')
% size(b1le)
% disp('b2le: ')
% size(b2le)
b3xi = b1le+b2le;
c3xi = c1le+c2le;

% //computing the trailing edge influence of the DVE element //

% //vortex system reference point midspan of DVE trailing edge
% //unless when neede for the alternative CDiEppler computation

if DVE_type == -4
    xsio(1) = 0.5*xsi;
else
    xsio(1) = xsi;
end

	a1te(1) = 0;  a1te(2) = 0; a1te(3) = 0;
	b1te(1) = 0;  b1te(2) = 0; b1te(3) = 0;
	c1te(1) = 0;  c1te(2) = 0; c1te(3) = 0;
    
    
if DVE_type == 0 || DVE_type == -2
    [a1te, b1te, c1te] = fcnBoundVortexInduction(Temp, xsiA, xsio, 0, phiTE, eta);
end

	a2te(1) = 0;  a2te(2) = 0; a2te(3) = 0;
	b2te(1) = 0;  b2te(2) = 0; b2te(3) = 0;
	c2te(1) = 0;  c2te(2) = 0; c2te(3) = 0;
    
if DVE_type ~= 3 && DVE_type ~= -3
    [a2te, b2te, c2te] = fcnVortexSheetInduction(Temp, xsiA, xsio, 0, phiTE, eta, singfct);
end

a3xi = a3xi - a1te;
b3xi = b3xi - (b1te + b2te);
c3xi = c3xi - (c1te + c2te);

% //rotating back into global co-system //
a3x = star_glob(a3xi,nu,eps,psi);
b3x = star_glob(b3xi,nu,eps,psi);
c3x = star_glob(c3xi,nu,eps,psi);

end

