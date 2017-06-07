function FW = fcnMove_Wing(FW,Temp,Aircraft)

%copy from freewake Move_Wing in wing_geometry.cpp

% //moves wing by delta x every time step,
% //function updates xo location of surface DVE's
% //
% //	input
% //	info			general information
% //	surfacePtr		surface DVE's
% //

%         //delta x = local U * delta time
delx(1) = Temp.u(1) * FW.Deltime;
delx(2) = Temp.u(2) * FW.Deltime;
delx(3) = Temp.u(3) * FW.Deltime;

for ii = 1:Aircraft.General.Panels % looping through influencing elements (all other DVEs)
    for jj = 1:length(FW.Panels(ii).DVE.Index)
        
        %         //move reference point
        FW.Panels(ii).DVE.xo(jj,:) = FW.Panels(ii).DVE.xo(jj,:) - delx;
        
        %move te points
        FW.Panels(ii).DVE.TECoordL(jj,:) = FW.Panels(ii).DVE.TECoordL(jj,:) -delx;
        FW.Panels(ii).DVE.TECoordR(jj,:) = FW.Panels(ii).DVE.TECoordR(jj,:) -delx;
    end
    %move TECoords
    FW.Panels(ii).TECoord(1,:) = FW.Panels(ii).TECoord(1,:) - delx; %left side
    FW.Panels(ii).TECoord(2,:) = FW.Panels(ii).TECoord(2,:) - delx; %right side
end