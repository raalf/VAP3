function [a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, xAi, xoj, nuj, phij, etaj, singfct)
% Facsimile of VortexSheetInduction from induced_velocity.cpp
% i is control point, j is inducing element

% ANGLES IN RADIANS!!! nuj, phij, etaj

%% From FreeWake comments:
% //this function computes influence coefficients of ctr point i due to
% //a semi-infinite vortex sheet.  The function is based on the function
% //FixedWakeInduction that is used in Horstmann's multiple lifting line
% //method.  The curren method has been modified to deal with the singularities
% //along the edges of the sheet, but is otherwise very similar to the function
% //FixedWakeInduction that is listed further below.
% //See also Appendix 3 of
% //"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
% //und  Nachrenchnung nichtplanarer Fluegelanordnungen",
% //by K.-H. Horstmann, 1987, DFVLR, Germany, DFVLR-FB 87-51.
% //
% //the return values are in global coordinates
% //
% //Input:
% //	xAi		element i control point in global coordinates
% //	x0j		element j midspan point of 1/4c line in global coordinates
% //	nuj		element j dihedral
% //	phij	element j sweep
% //	etaj	element j half span
% //  singfct rate at which singularity at edge of vortex sheet decays
% //			## singfct added 2/8/05 G.B.  This factor is only needed when
% //			relaxing the wake, otherwise it is set to zero (0)
% //
% //Output:
% //	a2x[3]	induced coefficient of bound vortex in global system
% //	b2x[3]	induced coefficient of bound vortex in global system
% //	c2x[3]	induced coefficient of bound vortex in global system

%% Function Body

xsij = abs(xoj(1));

if (singfct*singfct > 1000)
    disp('Issue with singfct in fcnVortexSheetInduction');
end

%   //transformation of i-control point global coordinates into
%   //local j-element coordinates. Horstmann Eq. 9
tempAA = -xoj;
tempA = tempAA + xAi;

if nuj*nuj > Temp.DBL_EPS
    tempAA = fcnRotateX(tempA,nuj);
    xsiA = tempAA(1);
    etaA = tempAA(2);
    zetaA = tempAA(3);
else
    xsiA = tempA(1);
    etaA = tempA(2);
    zetaA = tempA(3);
end

zetaASQR = zetaA*zetaA;

if zetaASQR <= Temp.DBL_EPS
    ABSzetaA = 0;
else
    ABSzetaA = sqrt(zetaASQR);
end

tempD = xsiA - etaA*tand(phij); % //if D=0 then A in plane of bound vortex
ABSD = abs(tempD); % // and zeta axis

if ABSzetaA <= Temp.DBL_EPS && ABSD <= Temp.DBL_EPS && abs(phij) <= Temp.DBL_EPS
    %//point A lies on bound vortex line/leading edge of vortexsheet
    
    b2x = [0 0 0]; c2x = [0 0 0];
    
    if abs(phij) <= Temp.DBL_EPS %//at leading edge of an unswept vortex sheet
        
        b2xi = [0 0]; c2xi = [0 0];
        
        t1 = etaA + etaj;
        t2 = etaA - etaj;
        
        % 		//Will be singular if t1 or t2 is 0,
        % 		//correction below requires that gamma=0 at tip!!!
        % 		//Thus, old version, changed 11/29/04
        % 		//old tempS = log((t1*t1)/(t2*t2));
        
        tempS = log((t1*t1 + singfct)/(t2*t2 + singfct));
        
        b2xi(3) = tempS*0.5;
        c2xi(3) = (-4*etaj + etaA*tempS);
        
        
        %//transform coefficients a, b, and c back into global co-system
        if nuj*nuj > Temp.DBL_EPS
            b2x = fcnRotateX(b2xi,-nuj);
            c2x = fcnRotateX(c2xi,-nuj);
        else
            b2x = b2xi;
            c2x = c2xi;
        end
    end
    
else
    % //point of interest does not fall on leading edge of sheet
    
    % //Horstmann Eq. A3-17
    a2 = 1 + tand(phij)*tand(phij);
    b2 = tempD*tand(phij);
    c2 = tempD*tempD + zetaASQR;
    
    t1 = etaA + etaj;
    t2 = etaA - etaj;
    
    rt1	= sqrt(t1*t1*a2 + 2*t1*b2 + c2);
    rt2	= sqrt(t2*t2*a2 + 2*t2*b2 + c2);
    
    %//Horstmann Eq. A3-18
    b(1) = -tempD;
    b(2) = zetaASQR*tand(phij);
    b(3) = 0;
    b(4) = -tand(phij);
    b(5) = -1;
    b(6) = 0;
    b(7) = 0;
    
    c(1) = -2*(b(2) - etaA*b(1));
    c(3) = 2*tand(phij);
    c(4) = 2*(xsiA - etaA*c(3));
    c(5) = -2*etaA;
    c(6) = -2*zetaASQR;
    c(7) = 2;
    c(2) = 0.5*c(6)*c(4);
    
    %//Horstmann Eq. A3-10
    epsilon = b(1)*b(1) + b(2)*b(4);
    rho = sqrt(epsilon*epsilon + 4*zetaASQR*b2*b2);
    
    tempS = 0.5*(rho + epsilon);
    
    if tempS <= Temp.DBL_EPS
        beta1 = 0;
    else
        beta1 = -sqrt(tempS);
    end
    
    tempS = 0.5*(rho-epsilon);
    
    if tempS <= Temp.DBL_EPS
        beta2 = 0;
    else
        beta2 = -sqrt(tempS);
    end
    
    %//see Horstmann program FLU lines 1443 through 1447
    tempS = zetaA*b2;
    if abs(tempS) > Temp.DBL_EPS
        beta2 = beta2*tempS/abs(tempS);
    end
    
    gamma1 	= (a2*beta2*zetaA + b2*beta1);
    gamma2 	= (a2*beta1*zetaA - b2*beta2);
    
    delta1 	= (b2*beta2*zetaA + c2*beta1);
    delta2 	= (b2*beta1*zetaA - c2*beta2);
    
    tempS	= gamma1*t1 + delta1 - rt1*rho;
    tempSS	= gamma2*t1 + delta2;
    
    mu1t1	= (tempS*tempS + tempSS*tempSS);
    
    %//mu2t1 computed as in Horstmann Eq. A3-10. atan-values then corrected.
    %atand down here????????, no, I checked (Bill).
    if t1*t1 <= Temp.DBL_EPS
        mu2t1 = pi/2*ABSzetaA/zetaA + atan(tempSS/tempS);
    else
        mu2t1 = atan(zetaA/t1) + atan(tempSS/tempS);
        %//Corrections according to Horstmamm program FLU lines 1473 and 1474
        if zetaA > 0 && t1 < 0
            mu2t1 = mu2t1 + pi;
        elseif zetaA < 0 && t1 < 0
            mu2t1 = mu2t1 - pi;
        end
    end
    
    %//Corrections according to Horstmamm program FLU lines 1484 and 1485
    if tempS < 0
        mu2t1 = mu2t1 + pi;
    elseif tempSS < 0 && tempS > 0
        mu2t1 = mu2t1 + 2*pi;
    end
    
    tempS = gamma1*t2 + delta1 - rt2*rho;
    tempSS = gamma2*t2 + delta2;
    
    mu1t2 = tempS*tempS + tempSS*tempSS;
    
    if t2*t2 < Temp.DBL_EPS
        mu2t2 = pi/2*ABSzetaA/zetaA + atan(tempSS/tempS); %// |zeta/t2| -> infinity
    else
        mu2t2 = atan(zetaA/t2) + atan(tempSS/tempS);
        % //Corrections according to Horstmamm program FLU lines 1466 1467
        if zetaA > 0 && t2 < 0
            mu2t2 = mu2t2 + pi;
        elseif zetaA < 0 && t2 < 0
            mu2t2 = mu2t2 - pi;
        end
    end
    
    %//Corrections according to Horstmamm program FLU lines 1479 and 1480
    if tempS < 0
        mu2t2 = mu2t2 + pi;
    elseif tempSS < 0 && tempS > 0
        mu2t2 = mu2t2 + 2*pi;
    end
    
    %     //Horstmann Eq. A3-13
    % 	//changed 11/23/04 G.B. in order to deal with leading edge singularity
    if abs(phij) <= Temp.DBL_EPS
        mu3t1 	= a2*t1 + b2 + sqrt(a2)*rt1;
        mu3t2 	= a2*t2 + b2 + sqrt(a2)*rt2;
    else
        mu3t1 	= .0001*xsij + a2*t1 + b2 + sqrt(a2)*rt1;
        mu3t2 	= .0001*xsij + a2*t2 + b2 + sqrt(a2)*rt2;
    end
    
    % 	//KHH eq. A3-8 through A3-16
    % 	//modified G25
    % 	//accounts for singularity along side edge of sheet
    % 	//used for b and c-influence coef. of zeta velocity
    % 	//added 4/13/04 G.B.
    % 	//added factor of 0.01*etaj in order to blend in faster singularity
    % 	//added 7/5/04 G.B.
    % 	//G25b includes the factor b25 (b[4])
    tempS	= singfct	+ zetaASQR + t1*t1;
    tempSS	= singfct	+ zetaASQR + t2*t2;
    
    G25b 	= -0.5*log(tempSS/tempS);
    
    % 	//G25c includes the factor c25 (c[4])		original - + - + (last one on line 1572)
    G25c 	= -etaj*log(tempS*tempSS);
    
    if abs(t1) > Temp.DBL_EPS %//point of interest not on left edge
        G25c = G25c + t1*log(zetaASQR + t1*t1);
    end
    
    if abs(t2) > Temp.DBL_EPS %//point of interest not on right edge
        G25c = G25c - t2*log(zetaASQR + t2*t2);
    end
    %     //point A in plane that is spaned by bound vortex and zeta-axis
    if ABSD <= Temp.DBL_EPS
        G2(1) = 0;
        G21b = 0;
        G21c = 0;
        G2(2) = 0;
        
        G2(5) = 0.5*log((t2*t2 + zetaASQR)/(t1*t1 + zetaASQR));
        %         		//G26 Horstmann Eq. A3-15, see also KHH program lines 1504 ff
        if ABSzetaA <= Temp.DBL_EPS
            G2(6) = 0;
        else %atand down here??????
            tempS = zetaASQR + t1*t2;
            G2(6) = atan((t2-t1)*zetaA/tempS);
            
            if tempS < 0 && (t2/zetaA) > 0
                G2(6) = G2(6) + pi;
            end
            
            if tempS <0 && (t2/zetaA) < 0
                G2(6) = G2(6) - pi;
            end
            
            G2(6) = G2(6)/zetaA;
        end
        
        
    else
        if ABSzetaA <= Temp.DBL_EPS
            % //point A lies i xsi-eta plane, but NOT on bound vortex line
            % 		 //another possible correction is required further down when
            % 	 	 //computing b2xi[1], c2xi[1],b2xi[2], c2xi[2]
            %
            % 	 	 //KHH line 1459 through 1492, removed 04/13/04 G.B.
            G2(1) = 0;
            %  			//modified G21, accounts for singularity along side edge of
            % 			//sheet used for b and c-influence coef. of zeta velocity
            % 			//include factors b21 and c21
            % 			//added 04/13/04 G.B. and modified 02/22.05
            tempS = log(mu1t2/mu1t1);
            G21b = b(1)*beta1*(0.5*tempS + G25b)/rho;
            G21c = b(1)*beta1*(etaA*tempS + G25c)/rho;
            
            G2(2) = 0;
            G2(6) = 0;
        else
            %             //point A is NEITHER in xsi-eta plane, NOR on bound vortex
            %
            % 			//G25 Horstmann Eq. A3-14
            % 			//used for b and c-influence coef. of eta velocity
            G2(5) = 0.5*log((t2*t2 + zetaASQR)/(t1*t1 + zetaASQR));
            
            tempS = 0.5*log(mu1t2/mu1t1);
            tempSS = mu2t2-mu2t1;
            %             //G21 Horstmann Eq. A3-8
            G2(1) = (beta1*(tempS-G2(5)) + beta2*tempSS)/rho;
            G21b = 0;
            G21c = 0;
            %             //G22 Horstmann Eq. A3-9
            G2(2) = (-beta2*(tempS-G2(5)) + beta1*tempSS)/(rho*zetaA);
            %             //G26 Horstmann Eq. A3-15, see also KHH program lines 1504 ff
            tempS = zetaASQR + t1*t2;
            G2(6) = atan((t2-t1)*zetaA/tempS);
            
            if tempS < 0 && (t2/zetaA) > 0
                G2(6) = G2(6) + pi;
            end
            
            if tempS <0 && (t2/zetaA) < 0
                G2(6) = G2(6) - pi;
            end
            
            G2(6) = G2(6)/zetaA;
        end
    end
    %     	//G24 Horstmann Eq. A3-12
    if abs(mu3t2) <= Temp.DBL_EPS || abs(mu3t1) <= Temp.DBL_EPS
        G2(4) = 0;
    else
        G2(4) = log(mu3t2/mu3t1)/sqrt(a2);
    end
    % 	//G23 Horstmann Eq. A3-11
    G2(3) = (rt2 - rt1 - b2*G2(4))/a2;
    % 	//G27 Horstmann Eq. A3-16
    G2(7) = t2-t1;
    
    % 	//Horstmann Eq. A3-7
    if ABSzetaA <= Temp.DBL_EPS
        %     //point A does lie in xsi-eta plane,
        % 		//wake can only induce zeta velocities
        %
        % 		//Horstmann Eq. A3-7, fist part b2eta, c2eta
        b2xi(2) = 0;
        c2xi(2) = 0;
        %    		//Horstmann Eq. A3-7, second part bzeta, czeta
        b2xi(3) = G21b + G2(4)*b(4) + G25b;
        c2xi(3) = G21c + G2(3)*c(3) + G2(4)*c(4) + G25c + G2(7)*c(7);
        
    else
        % //point does NOT lie in xsi-eta plane, zetaA NOT= 0
        % 		//THIS IS THE NON-EXCEPTION
        % //Horstmann Eq. A3-7, fist part b2eta, c2eta
        b2xi(2) = -zetaA*(G2(1)*b(4) + G2(2)*b(1) + G2(6)*b(5));
        c2xi(2) = -zetaA*(G2(1)*c(4) + G2(2)*c(1) + G2(4)*c(3) + G2(5)*c(7) + G2(6)*c(5));
        %     //Horstmann Eq. A3-7, second part bzeta, czeta
        b2xi(3) = G2(1)*b(1) + G2(2)*b(2) + G2(4)*b(4) + G2(5)*b(5);
        c2xi(3) = G2(1)*c(1) + G2(2)*c(2) + G2(3)*c(3) + G2(4)*c(4) + G2(5)*c(5) + G2(6)*c(6) + G2(7)*c(7);
        
    end
    
    %     //wake only induces velocities in eta and zeta direction
    %     //CAUTION: that requires small angles or xsi and Uinf aligned
    b2xi(1) = 0;
    c2xi(1) = 0;
    
    if nuj*nuj > Temp.DBL_EPS
        
        b2x = fcnRotateX(b2xi,-nuj);
        c2x = fcnRotateX(c2xi,-nuj);
    else
        b2x = b2xi;
        c2x = c2xi;
        
    end

    
end
    % //a2 coefficients of vortex sheet are zero, since vorticity is B+2*C*eta
    a2x = [0 0 0];
end

