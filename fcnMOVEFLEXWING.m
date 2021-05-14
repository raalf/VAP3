function [SURF, MISC, COND] = fcnMOVEFLEXWING(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, valTIMESTEP)
% MISC.matNEWWAKE, MISC.matNPNEWWAKE,
% This function determines the velocities with which the DVEs are moved
% based on the deflection and twist of the wing. The corresponding
% translations are then computed of the DVE vertices and control points.

% if COND.valFULLTRIMSTEP > 1 && FLAG.FULLTRIM == 0
%     
%     OUTP.matDEFGLOB = OUTP.matDEFGLOBTRIM;
%     OUTP.matTWISTGLOB = OUTP.matTWISTGLOBTRIM;
%     
% end

[ledves, ~, ~] = find(SURF.vecDVELE(SURF.idxFLEX) > 0);

% Span of each spanwise set of DVEs
vecDVESPAN = 2*SURF.vecDVEHVSPN(ledves)';

% Calculate cartesian velocity of DVE edges

matUINF_edge = [interp1(SURF.matCENTER(ledves,2),SURF.matUINF(ledves,1),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.matCENTER(ledves,2),SURF.matUINF(ledves,2),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.matCENTER(ledves,2),SURF.matUINF(ledves,3),SURF.vecSPANLOC,'linear','extrap')];

del_twist = ((OUTP.matTWISTGLOB(valTIMESTEP,:) - OUTP.matTWISTGLOB(valTIMESTEP-1,:)));
omega = ((OUTP.matTWISTGLOB(valTIMESTEP,:) - OUTP.matTWISTGLOB(valTIMESTEP-1,:)))./COND.valDELTIME;
vecXVEL = matUINF_edge(:,1);
vecYVEL = matUINF_edge(:,2) + [0, (OUTP.matDEFGLOB(valTIMESTEP,2:end)-OUTP.matDEFGLOB(valTIMESTEP-1,2:end))./...
    (COND.valDELTIME.*tan(repmat(pi/2,1,size(OUTP.matSLOPE,2)-1)-(OUTP.matSLOPE(valTIMESTEP,2:end)-OUTP.matSLOPE(valTIMESTEP-1,2:end))./2))]';
vecZVEL = matUINF_edge(:,3) + ((OUTP.matDEFGLOB(valTIMESTEP,:) - OUTP.matDEFGLOB(valTIMESTEP-1,:))./COND.valDELTIME)';

% Determine DVEs in each spanwise station
[matROWS] = fcnDVEROW(ledves, SURF, INPU, 1);

%% Determining displacement distances

% Determine vertices that need to be moved at each spanwise station

% All left LE and TE points to move
temp_leftV = [SURF.matNPDVE(matROWS{1},1),SURF.matNPDVE(matROWS{1},4)];
temp_leftV = reshape(temp_leftV,sum(INPU.vecN(FLAG.vecFLEXIBLE == 1),1),[]);

[move_row,~] = find(temp_leftV); % Vector correspond to which index of deflection velocity matrix should be used for each element

tempENDWING = max(max(SURF.matNPDVE(SURF.idxFLEX,:))); % Index for last row in matNPVLST that corresponds to the wing
tempENDTAIL = max(max(SURF.matNPDVE(SURF.idxTAIL,:)));

% Allocate space for translation matrices
translateNPVLST = zeros(size(SURF.matNPVLST(1:tempENDWING),1),3);
temp_translate = zeros(size(SURF.matNPVLST(1:tempENDWING),1),3);

temp_r = [sqrt(sum(SURF.matSCLST.^2,2)), zeros(length(SURF.matSCLST),2)]; % Distance between vertex and shear center

xz_sign = sign(SURF.matSCLST(1:tempENDWING,1)); % Determines whether positive or negative contribution to X and Z velocity for each vertex

% Perform a linear interpolation/extrapolation to determine the pitch of
% the DVE's at their respective left and right edges. This is used to
% determine the initial orientation of the DVE before applying the twist
% caused by elastic deformation
tempSPANDIST = SURF.matCENTER(ledves,2); % Y coordinate of DVE mid-point (point where DVEPITCH is applied) --> used as "x" term for linear interpolation

vecEDGEPITCH = SURF.vecDVEPITCH(ledves); % DVE pitch along span --> used as "y" term for linear interpolation

% Some setup of work to be able to perform linear interpolation without a
% for loop
tempSPANDIST = repmat(tempSPANDIST', size(tempSPANDIST,1),1);

tempSPANDIST = triu(tempSPANDIST);

vecEDGEPITCH = repmat(vecEDGEPITCH', size(vecEDGEPITCH,1),1);

vecEDGEPITCH = triu(vecEDGEPITCH);

vecEDGEPITCH = ((SURF.vecSPANLOC(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
    tempSPANDIST(1,1:(end-1)))).*(vecEDGEPITCH(2,2:end)-vecEDGEPITCH(1,1:(end-1))) + vecEDGEPITCH(1,1:(end-1)); % Linear interpolation

% Adding in root and tip values using a linear extrapolation
pitch_root = vecEDGEPITCH(1,1) - (tempSPANDIST(1,1) - SURF.vecSPANLOC(1)).*(vecEDGEPITCH(1,2)-vecEDGEPITCH(1,1))./(tempSPANDIST(1,2)-tempSPANDIST(1,1));

pitch_tip = vecEDGEPITCH(1,end) + (SURF.vecSPANLOC(end) - tempSPANDIST(1,end)).*(vecEDGEPITCH(1,end)-vecEDGEPITCH(1,end-1))./(tempSPANDIST(1,end)-tempSPANDIST(1,end-1));

vecEDGEPITCH = [pitch_root, vecEDGEPITCH, pitch_tip];

% ======================== Left Edge Displacements ========================
% Translate left edge vertices due to twist
twistXDIST = -temp_r(temp_leftV,1).*cos(del_twist(move_row)+vecEDGEPITCH(move_row))' + ...
    temp_r(temp_leftV,1).*cos(vecEDGEPITCH(move_row))'; % X component of twist 
twistZDIST = temp_r(temp_leftV,1).*sin(del_twist(move_row)+vecEDGEPITCH(move_row))' -  ...
    temp_r(temp_leftV,1).*sin(vecEDGEPITCH(move_row))'; % Z component of twist

v_rot = temp_r(temp_leftV,1).*omega(move_row)';

% Assign twist displacement to translation matrix
temp_translate(temp_leftV,1) = twistXDIST;
temp_translate(temp_leftV,3) = twistZDIST;

test = zeros(size(SURF.matNPVLST,1),3);
test(temp_leftV,1) = v_rot.*sin(OUTP.matTWISTGLOB(valTIMESTEP,move_row)+vecEDGEPITCH(move_row))'.*COND.valDELTIME;
test(temp_leftV,3) = v_rot.*cos(OUTP.matTWISTGLOB(valTIMESTEP,move_row)+vecEDGEPITCH(move_row))'.*COND.valDELTIME;

% Translate left edge vertices due to freestream and bending
translateNPVLST(temp_leftV,1) = COND.valDELTIME.*vecXVEL(move_row);
translateNPVLST(temp_leftV,2) = COND.valDELTIME.*vecYVEL(move_row);
translateNPVLST(temp_leftV,3) = -1*COND.valDELTIME.*vecZVEL(move_row);

% ======================== Right Edge Displacements =======================
% All right LE and TE points to move
temp_rightV = [SURF.matNPDVE(matROWS{1},2), SURF.matNPDVE(matROWS{1},3)];
temp_rightV = reshape(temp_rightV,sum(INPU.vecN(FLAG.vecFLEXIBLE == 1),1),[]);

[move_row,~] = find(temp_rightV); % Vector correspond to which index of deflection velocity matrix should be used for each element

v_rot = temp_r(temp_leftV,1).*omega(move_row+1)';
% Translate right edge vertices due to twist
twistXDIST = -temp_r(temp_rightV,1).*cos(del_twist(move_row+1)+vecEDGEPITCH(move_row+1))' + ...
    temp_r(temp_rightV,1).*cos(vecEDGEPITCH(move_row+1))'; % X component of twist
twistZDIST = temp_r(temp_rightV,1).*sin(del_twist(move_row+1)+vecEDGEPITCH(move_row+1))' - ...
    temp_r(temp_rightV,1).*sin(vecEDGEPITCH(move_row+1))'; % Z component of twist 

temp_translate(temp_rightV,1) = twistXDIST';
temp_translate(temp_rightV,3) = twistZDIST';

test(temp_rightV,1) = v_rot.*sin(OUTP.matTWISTGLOB(valTIMESTEP,move_row+1)+vecEDGEPITCH(move_row+1))'.*COND.valDELTIME;
test(temp_rightV,3) = v_rot.*cos(OUTP.matTWISTGLOB(valTIMESTEP,move_row+1)+vecEDGEPITCH(move_row+1))'.*COND.valDELTIME;

% Assign appropriate sign to twist movement
temp_translate(:,1) = xz_sign.*temp_translate(:,1);
temp_translate(:,3) = xz_sign.*temp_translate(:,3);

% ======================================================================= %
% =============== TRYING NEW TRANSLATION FROM TWIST ===================== %
% ======================================================================= %
% temp_translate(:,1) = xz_sign.*test(:,1);
% temp_translate(:,3) = xz_sign.*test(:,3);
% ======================================================================= %
% ======================================================================= %

% Translate right edge vertices due to freestream and bending
translateNPVLST(temp_rightV,1) = COND.valDELTIME.*vecXVEL(move_row+1);
translateNPVLST(temp_rightV,2) = COND.valDELTIME.*vecYVEL(move_row+1);
translateNPVLST(temp_rightV,3) = -1*COND.valDELTIME.*vecZVEL(move_row+1);

if any(FLAG.vecTRIMABLE == 1) == 1
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matVEHUVW*COND.valDELTIME;
    
end

%% Move wing and generate new wake elements

% Old trailing edge vertices
MISC.matNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,4) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,3) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),4),:);
MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),3),:);

% Update SURF.matVLST and SURF.matNTVLST
SURF.matNPVLST(1:tempENDWING,:) = SURF.matNPVLST(1:tempENDWING,:) - (translateNPVLST - temp_translate);

% Move stiff tail if it exists
SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matVEHUVW,size(SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:),1),1);
SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matVEHUVW,size(SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:),1),1);

% New non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

