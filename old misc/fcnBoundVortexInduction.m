function [a1x, b1x, c1x] = fcnBoundVortexInduction(Temp, xAi, xoj, nuj, phij, etaj)
% Facsimile of function "BoundVortexInduction" from induced_velocity.cpp
% i is the current element, j is the influencing element

% ANGLES IN RADIANS!!! nuj, phij, etaj

%% From FreeWake comments:
% //this function computes influence coefficients of ctr point i due to
% //the bound vortex at elementary wing j.
% //See also Appendix 3 of
% //"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
% //und  Nachrenchnung nichtplanarer Fluegelanordnungen",
% //by K.-H. Horstmann, published 1987 by DFVLR, Germany, DFVLR-FB 87-51.
%
% //the return values are in global coordinates
% //
% //Input:
% //	xAi		control point i in global coordinates
% //	x0j		element j midspan point of 1/4c line in global coordinates
% //	nuj		element j dihedral
% //	phij	element j sweep
% //	etaj	element j half span
% //
% //Output:
% //	a1x[3]	induced coefficient of bound vortex in global system
% //	b1x[3]	induced coefficient of bound vortex in global system
% //	c1x[3]	induced coefficient of bound vortex in global system

%% Function Body

tempAA = -xoj;
tempA = tempAA + xAi;

if nuj*nuj > Temp.DBL_EPS
    clear tempAA
    tempAA = fcnRotateX(TempA, nuj);
    
    xsiA = tempAA(1);
    etaA = tempAA(2);
    zetaA = tempAA(3);
    
else
    xsiA = tempA(1);
    etaA = tempA(2);
    zetaA = tempA(3);
end

% 	//if D=0 then A in plane of bound vortex and zeta axis
tempD = xsiA-etaA*tand(phij);

if abs(zetaA) < Temp.ZERO && abs(tempD) < Temp.ZERO
    
    %       //1. condition: point A lies in xsiA-etaA-plane
    % 		//	=>only velocities in zeta-direction induced
    % 		//2. condition: point A lies in plane defined by zeta-axis and
    % 		// 	 bound vortex line
    % 		//  =>no velocity in zeta-direction induced
    % 		//or in other words, point lies on line of vortex filament
    
    a1xi = [0 0 0]; b1xi = [0 0 0]; c1xi = [0 0 0];
    
else
    
    a1 = 1 + tand(phij)*tand(phij);
    b1 = -(etaA + xsiA*tand(phij));
    c1 = xsiA*xsiA + etaA*etaA + zetaA*zetaA;
    
    reta1 = sqrt(etaj*etaj*a1 - 2*etaj*b1 + c1);
    reta2 = sqrt(etaj*etaj*a1 + 2*etaj*b1 + c1);
    
    tempS = ((a1*c1 - b1*b1)*reta1*reta2);
    G11	 = 	(a1*etaj*(reta1 + reta2) + b1*(reta1 - reta2))/tempS;
    G12	 = 	(-b1*etaj*(reta1 + reta2) - c1*(reta1 - reta2))/tempS;
    G13	 = 	((2*b1*b1 - a1*c1)*etaj*(reta1 + reta2) + b1*c1*(reta1 - reta2))/(a1*tempS);
    G13 = 	G13 + log((sqrt(a1)*reta2 + a1*etaj + b1)/(sqrt(a1)*reta1 - a1*etaj + b1))/sqrt(a1*a1*a1);
    
    % // Horstmann Eq. A3-2
    if abs(zetaA) < Temp.ZERO
        %       //condition: point A lies in xsiA-etaA-plane
        % 		//	=>only velocities in zeta-direction induced
        a1xi = [0 0]; b1xi = [0 0]; c1xi = [0 0];
        
    else
        a1xi(1) = -G11*zetaA;
        a1xi(2) = -a1xi(1)*tand(phij);
        
        b1xi(1) = a1xi(1)*G12/G11;
        b1xi(2) = a1xi(2)*G12/G11;
        
        c1xi(1) = a1xi(1)*G13/G11;
        c1xi(2) = a1xi(2)*G13/G11;
        
    end
    
    if abs(tempD) < Temp.ZERO
        
        %       //condition: point A lies in plane defined by zeta-axis and
        %       //bound vortex line
        %       // =>no velocity in zeta-direction induced
        a1xi(3) = 0; b1xi(3) = 0; c1xi(3) = 0;
        
    else
        a1xi(3) = G11*tempD;
        b1xi(3) = a1xi(3)*G12/G11;
        c1xi(3) = a1xi(3)*G13/G11;
        
    end
end

% //transform coefficients a, b, and c back into global co-system
if nuj*nuj > Temp.DBL_EPS
    a1x = fcnRotateX(a1xi,-nuj);
    b1x = fcnRotateX(b1xi,-nuj);
    c1x = fcnRotateX(c1xi,-nuj);
else
    a1x = a1xi; b1x = b1xi; c1x = c1xi;
    
end



end

