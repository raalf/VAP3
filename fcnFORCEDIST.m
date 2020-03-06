function [OUTP] = fcnFORCEDIST(SURF, COND, INPU, OUTP, FLAG, valTIMESTEP)

% This function computes the dimensional force and moment distribution
% across the wing, resolved to the shear center line. Moment is taken
% about the elastic axis.
%
% INPUT:
%
% OUTPUT:
% OUTP.vecLIFTDIST - 1 x sum(INPU.vecN) matrix of the total lift at each spanwise
% station
% OUTP.vecMOMDIST - 1 x sum(INPU.vecN) matrix of the total pitching moment at each
% spanwise station (+ve nose up)

% Calculate qinf required for steady level flight. This will be used for
% load calculations
% q_inf = valWEIGHT/(valCL*valAREA);
% valVINF = sqrt(2*q_inf/COND.valDENSITY);
q_inf = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF;

[ledves, ~, ~] = find(SURF.vecDVELE(SURF.idxFLEX) > 0);

lepanels = SURF.vecDVEPANEL(ledves);

% Determine DVEs in each spanwise station
for i = 1:length(find(FLAG.vecFLEXIBLE == 1))

	idxdve = ledves(SURF.vecDVEWING(ledves) == i);
	idxpanel = lepanels(SURF.vecDVEWING(ledves) == i);

    m = INPU.vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
    
    matROWS{i} = repmat(idxdve,1,m) + double(tempm);

end

% ======================================================================= % 
%                        Start Lift Calculation                           %
% ======================================================================= %

% Convert DVE force from force/density to force/unit length
% OUTP.vecLIFTDIST = ((sum(SURF.vecDVELFREE(matROWS{1}),2) + sum(SURF.vecDVELIND(matROWS{1}),2))*COND.valDENSITY)./(2*SURF.vecDVEHVSPN(ledves));
% 
% % Interpolate lift distribution to find lift at DVE edges (i.e. Structural
% % grid points)
% tempSPANDIST = SURF.matCENTER(ledves,2); % Y coordinate of DVE midpoints
% 
% % Some setup work to be able to perform linear interpolation without a
% % for loop
% tempSPANDIST = repmat(tempSPANDIST', size(tempSPANDIST,1),1);
% 
% tempSPANDIST = triu(tempSPANDIST);
% 
% OUTP.vecLIFTDIST = repmat(OUTP.vecLIFTDIST', size(OUTP.vecLIFTDIST,1),1);
% 
% OUTP.vecLIFTDIST = triu(OUTP.vecLIFTDIST);
% 
% OUTP.vecLIFTDIST = ((SURF.vecSPANDIST(2:(end-1)) - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
%     tempSPANDIST(1,1:(end-1)))).*(OUTP.vecLIFTDIST(2,2:end)-OUTP.vecLIFTDIST(1,1:(end-1))) + OUTP.vecLIFTDIST(1,1:(end-1)); % Linear interpolation of lift
% 
% OUTP.vecLIFTDIST = [OUTP.vecLIFTDIST(1), OUTP.vecLIFTDIST, 0]'; % Add zero lift to tip. Use the same lift at root as at the next neighbouring node. This doesn't really have any impact since the structure solver doesn't use this value
% 
% OUTP.vecWRBM(valTIMESTEP,1) = sum(SURF.vecSPANDIST.*OUTP.vecLIFTDIST,1);

% ======================================================================= % 
%                        Start Moment Calculation                         %
% ======================================================================= %
lift_chord = SURF.vecDVELFREE(matROWS{1}) + SURF.vecDVELIND(matROWS{1}); % Lift distribution along each chordwise location

% Interpolate lift force at each chordwise station to use for moment
% calculations. Also determine the LE coordinates of each chordwise DVE to
% be used to calculate a moment arm between DVE and elastic axis
for i = 1:INPU.vecM(1)
    
    row_ledves(:,:,i) = [SURF.matNPVLST(SURF.matNPDVE(matROWS{1}(:,i),1),:); SURF.matNPVLST(SURF.matNPDVE(matROWS{1}(end,i),2),:)]; % LE coordinates of each set of chordwise DVEs
    
    temp_liftdist = repmat(lift_chord(:,i)', size(lift_chord,1),1);

    temp_liftdist = triu(temp_liftdist);

    temp_liftdist = ((SURF.vecSPANDIST(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
        tempSPANDIST(1,1:(end-1)))).*(temp_liftdist(2,2:end)-temp_liftdist(1,1:(end-1))) + temp_liftdist(1,1:(end-1)); % Linear interpolation of lift
    
    temp_lift = temp_liftdist;
   
    temp_lift = [temp_lift(1); temp_liftdist']; % Use the same lift at root as at the next neighbouring node
    
    % Lift force in X, Y, Z components to be used in cross product to
    % calculate moment about elastic axis. Each step into the 3rd dimension
    % is the lift distribution at each chordwise station
    lift_moment(:,:,i) = [temp_lift.*SURF.matLIFTDIR(matROWS{1}(:,i),:); zeros(1,3)]; 
    
end

SURF.matSC = repmat(SURF.matSC,1,1,INPU.vecM(1));

% Compute moment arm for cross product
% delX = (row_ledves - SURF.matSC);
delX = SURF.matAEROCNTR - SURF.matSC;
% delX = [vecLSAC, zeros(size(vecLSAC,1),2)];
% force = [zeros(size(OUTP.vecLIFTDIST,2),2),OUTP.vecLIFTDIST'];
tempMOMDIST = cross(delX,lift_moment); % M' = delx X lift
% tempMOMDIST = OUTP.vecLIFTDIST'.*vecLSAC;

% Compute magnitude of moment at each chordwise location
for i = 1:INPU.vecM(1)
   
    matMOMDIST(:,i) = sign(tempMOMDIST(:,2,i)).*sqrt(sum(abs(tempMOMDIST(:,:,i)).^2,2));
    
end

% vecGAMMA = A(SURF.vecLEDVES) + SURF.vecDVEHVSPN(SURF.vecLEDVES)'.*B(SURF.vecLEDVES) + SURF.vecDVEHVSPN(SURF.vecLEDVES)'.*SURF.vecDVEHVSPN(SURF.vecLEDVES)'.*C(SURF.vecLEDVES);
% gamma_root = A(SURF.vecLEDVES(1)) - SURF.vecDVEHVSPN(SURF.vecLEDVES(1))'.*B(SURF.vecLEDVES(1)) + SURF.vecDVEHVSPN(SURF.vecLEDVES(1))'.*SURF.vecDVEHVSPN(SURF.vecLEDVES(1))'.*C(SURF.vecLEDVES(1));
% 
% vecGAMMA = [gamma_root,vecGAMMA];

OUTP.vecMOMDIST = sum(matMOMDIST,2);
OUTP.vecMOMDIST = [OUTP.vecMOMDIST(1:(end-1)).*COND.valDENSITY./(2*SURF.vecDVEHVSPN(SURF.vecLEDVES)); 0];
% M0 = q_inf*valAREA.*SURF.vecMAC'.*SURF.vecMAC'.*valCM./vecGAMMA(1:(end-1));

% M = [s(SURF.vecLEDVESm,:).*M0'; zeros(1,3)];

% OUTP.vecMOMDIST = OUTP.vecMOMDIST + M(:,2);

OUTP.vecMOMDIST = OUTP.vecMOMDIST + [q_inf*SURF.vecMAC.*SURF.vecMAC.*OUTP.vecCMDIST;0];

end