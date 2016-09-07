function [w_ind] = fcnSingle_DVE_Induced_Velocity(FW,Temp,timestep,Panel,DVE,P,DVE_type,Surf_type)

% Facsimile of Single_DVE_Induced_Velocity from induced_velocity.cpp

% Panel is the Panel ##, DVE is the DVE ## in that panel
%Surf_type is 1 for Surface, 2 for wake
%% From FW Comments:
% //computes induced velocity in point P due to DVE DVelement.
% //also computes induced velocity of symmetry element
% //
% //input:
% //	info		general info
% //	DVelement	distributed vorticity element that induces
% //	P			point of interest
% // 	DVE_type 	type of DVE,
% //	 	DVE_type = 0 DVE has vortex filaments at leading and
% //		 		 	trailing edge, usually lifting surface DVE
% //		DVE_type = 1 DVE has no vortex filaments at leading and
% //					trailing edge, usually wake DVE
% //		DVE_type = 2 DVE has vortex filament at leading edge,
% //					but not at trailing edge
% //		DVE_type =-2 DVE has vortex filament at trailing edge,
% //					but not at leading edge
% //		DVE_type = 3 DVE is a semi infinite vortex sheet without a
% //					vortex filaments at leading and trailing edge
% //		DVE_type =-3 DVE is a semi infinite vortex sheet with a
% //					 vortex filaments at its leading edge
% //		DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
% //					 xsi aft of the ref. pt. (for CD computation along TE)
% //		DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
% //					 also vortex filament at LE (CD computation along TE)
%
% //
% //
% //output:
% //	w_ind		induced velocity in P due to DVelement

% 	//computes influence coefficients, a3,b3, and c3, of
% 	//DVE on on point P.
% 	//(symmetry case is considered in subroutine)


%for the wake, timestep + 1 is the current timestep's index, since timestep
%0 is the first row of wake elements. 


if Surf_type == 1 %surface DVE
    singfct = FW.Panels(Panel).DVE.singfct(DVE); %CHANGE TO SURFACE ELEMENTS!!! % Changed 2016-02-07, T.D.K
else
    %Surf_type == 2 %wake DVE
    singfct = FW.Panels(Panel).WakeDVE(timestep + 1).singfct(DVE);
end

[a3, b3, c3] = fcnDVEInfluenceCoeff(FW, Temp, timestep, FW.Sym, Panel, DVE, P, singfct, DVE_type, Surf_type);

%     	//velocities (sort of like KHH eq. 36 and 38)
if Surf_type == 1 %surface DVE
    w_ind(1) = FW.Panels(Panel).DVE.A(DVE)*a3(1) + FW.Panels(Panel).DVE.B(DVE)*b3(1) + FW.Panels(Panel).DVE.C(DVE)*c3(1);
    w_ind(2) = FW.Panels(Panel).DVE.A(DVE)*a3(2) + FW.Panels(Panel).DVE.B(DVE)*b3(2) + FW.Panels(Panel).DVE.C(DVE)*c3(2);
    w_ind(3) = FW.Panels(Panel).DVE.A(DVE)*a3(3) + FW.Panels(Panel).DVE.B(DVE)*b3(3) + FW.Panels(Panel).DVE.C(DVE)*c3(3);
else
    %Surf_type == 2 %wake DVE
    w_ind(1) = FW.Panels(Panel).WakeDVE(timestep + 1).A(DVE)*a3(1) + FW.Panels(Panel).WakeDVE(timestep + 1).B(DVE)*b3(1) + FW.Panels(Panel).WakeDVE(timestep + 1).C(DVE)*c3(1);
    w_ind(2) = FW.Panels(Panel).WakeDVE(timestep + 1).A(DVE)*a3(2) + FW.Panels(Panel).WakeDVE(timestep + 1).B(DVE)*b3(2) + FW.Panels(Panel).WakeDVE(timestep + 1).C(DVE)*c3(2);
    w_ind(3) = FW.Panels(Panel).WakeDVE(timestep + 1).A(DVE)*a3(3) + FW.Panels(Panel).WakeDVE(timestep + 1).B(DVE)*b3(3) + FW.Panels(Panel).WakeDVE(timestep + 1).C(DVE)*c3(3);
end
w_ind = w_ind*-1/(4*pi);


end