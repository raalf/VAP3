function [a, b, c] = fcnDVEInfluenceCoeff(FW, Temp, timestep, sym_condition, Panel, DVE, P, singfct, DVE_type, Surf_type)

% Facsimile of DVE_Influence_Coeff from induced_velocity.cpp

% Panel is the Panel ##, DVE is the DVE ## in that panel

% WHAT ABOUT SYMMETRY???

%% From FW Comments:

% //computes the influence coefficient due to the induction of a
% //distributed vortex element at point P. Considers symmetrical conditions.
% //input:
% 	// element	DV element that induces
% 	// P		point in which DVE induces
% 	// info		general setting info, such as symmetry
% 	// DVE_type 	type of DVE,
% 	//	 DVE_type = 0 DVE has vortex filaments at leading and
% 	//			 	trailing edge, usually lifting surface DVE
% 	//	 DVE_type = 1 DVE has no vortex filaments at leading and
% 	//				trailing edge, usually wake DVE
% 	//	 DVE_type = 2 DVE has vortex filament at leading edge,
% 	//				but not at trailing edge
% 	//	 DVE_type =-2 DVE has vortex filament at trailing edge,
% 	//				but not at leading edge
% 	//	 DVE_type = 3 DVE is a semi infinite vortex sheet without a
% 	//				 vortex filaments at leading and trailing edge
% 	//	DVE_type =-3 DVE is a semi infinite vortex sheet with a
% 	//				 vortex filaments at its leading edge
% 	//	DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
% 	//			 xsi aft of the ref. pt. (for CD computation along TE)
% 	//	DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
% 	//				also vortex filament at LE (CD computation along TE)
% 	//
% 	//
% 	//output:
% 	//a, b, c	influence coefficients as described by KHH in Appendix 3

%% Function Body:
if Surf_type == 1 %surface DVE
    xo = FW.Panels(Panel).DVE.xo(DVE,:);
    nu = FW.Panels(Panel).DVE.roll(DVE);
    eps = FW.Panels(Panel).DVE.pitch(DVE);
    phiLE = FW.Panels(Panel).DVE.phiLE(DVE);
    phiTE = FW.Panels(Panel).DVE.phiTE(DVE);
    psi = FW.Panels(Panel).DVE.yaw(DVE);
    eta = FW.Panels(Panel).DVE.eta(DVE);
    xsi = FW.Panels(Panel).DVE.xsi(DVE);
   %singfct is passed into this function since DVEKinCond needs to send 0
    
else
    %Surf_type == 2 %wake DVE
    xo = FW.Panels(Panel).WakeDVE(timestep+1).xo(DVE,:);
    nu = FW.Panels(Panel).WakeDVE(timestep+1).roll(DVE);
    eps = FW.Panels(Panel).WakeDVE(timestep+1).pitch(DVE);
    phiLE = FW.Panels(Panel).WakeDVE(timestep+1).phiLE(DVE);
    phiTE = FW.Panels(Panel).WakeDVE(timestep+1).phiTE(DVE);
    psi = FW.Panels(Panel).WakeDVE(timestep+1).yaw(DVE);
    eta = FW.Panels(Panel).WakeDVE(timestep+1).eta(DVE);
    xsi = FW.Panels(Panel).WakeDVE(timestep+1).xsi(DVE);
     %singfct is passed into this function since DVEKinCond needs to send 0
    
end
%fprintf('P = %f\t%f\t%f\n',P(1),P(2),P(3));

%nu, eps, phile and phite go in as deg
[a3, b3, c3] = fcnDVEInduction(Temp, P, xo, nu, eps, phiLE, phiTE, psi, eta, xsi, DVE_type, singfct);

a = a3; b = b3; c = c3;

% fprintf('a3: %f %f %f\n',a(1),a(2),a(3));
% fprintf('b3: %f %f %f\n',b(1),b(2),b(3));
% fprintf('c3: %f %f %f\n',c(1),c(2),c(3));

if sym_condition == 1
    tempA = [xo(1) -xo(2) xo(3)];
    
    [a3, b3, c3] = fcnDVEInduction(Temp, P, tempA, -nu, eps, -phiLE, -phiTE, -psi, eta, xsi, DVE_type, singfct);
    %     disp(b3)
    a = a + a3;
    b = b - b3;
    c = c + c3;
    
end
end

