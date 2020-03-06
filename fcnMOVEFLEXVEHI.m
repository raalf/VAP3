function [SURF, MISC, COND, INPU, VEHI] = fcnMOVEFLEXVEHI(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP)

tempENDWING = max(max(SURF.matNPDVE(SURF.idxFLEX,:))); % Index for last row in matNTVLST that corresponds to the wing
tempENDTAIL = max(max(SURF.matNPDVE(SURF.idxTAIL,:)));

% Deflection and twist at the wing control point y-locations
cpt_def(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP,:),SURF.center_dist);
cpt_twist(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP,:),SURF.center_dist);

% Compute deflection at propeller spanwise location if it exists
if isempty(VEHI.vecPROPLOC) == 0
    prop_def(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP,:),VEHI.vecPROPLOC_START(:,2));
    prop_twist(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP,:),VEHI.vecPROPLOC_START(:,2));
    prop_def(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP-1,:),VEHI.vecPROPLOC_START(:,2));
    prop_twist(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP-1,:),VEHI.vecPROPLOC_START(:,2));
end

cpt_def(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP-1,:),SURF.center_dist);
cpt_twist(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP-1,:),SURF.center_dist);

tempEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),SURF.vecWINGCG(:,2));
tempEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.vecWINGCG(:,2));
tempEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),SURF.vecWINGCG(:,2));

% Translation matrix for NTVLST points
for i = 1:size(SURF.idx_struct,2)

    % Perform rotation due to wing twist about elastic axis for DVE pts
    tempNPVLST = SURF.matNPVLST(SURF.idx_struct(:,i),:) - SURF.matEALST(SURF.idx_struct(:,i),:);

    deltaEPS = (OUTP.matTWISTGLOB(valTIMESTEP,i)-OUTP.matTWISTGLOB(valTIMESTEP-1,i));
    ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];
    
    vlst2 = (ROT*tempNPVLST')' + SURF.matEALST(SURF.idx_struct(:,i),:);
    
    rotNPVLST(SURF.idx_struct(:,i),:) = vlst2 - SURF.matNPVLST(SURF.idx_struct(:,i),:);
    
    % Perform rotation due to wing twist about elastic axis for wing CG
    if i > 1
        tempCG = (SURF.vecWINGCG(i-1,:) - tempEA(i-1,:))*0;

        deltaEPSCG = (cpt_twist(valTIMESTEP,i-1)-cpt_twist(valTIMESTEP-1,i-1));
        ROTCG = [cos(deltaEPSCG) 0 sin(deltaEPSCG); 0 1 0; -sin(deltaEPSCG) 0 cos(deltaEPSCG)];

        cg2 = (ROTCG*tempCG')' + SURF.vecWINGCG(i-1,:);

        rotCG(i-1,:) = cg2 - SURF.vecWINGCG(i-1,:);
        
        SURF.vecWINGCG(i-1,3) = SURF.vecWINGCG(i-1,3) + VEHI.matVEHUVW(:,3).*COND.valDELTIME + (cpt_def(valTIMESTEP,i-1) - cpt_def(valTIMESTEP-1,i-1)) + rotCG(i-1,3);
        SURF.vecWINGCG(i-1,1) = SURF.vecWINGCG(i-1,1) + VEHI.matVEHUVW(:,1).*COND.valDELTIME + rotCG(i-1,1);
    end
    
    % Move DVE pts and elastic axis in x and z direction based on wing
    % bending and twist
    SURF.matNPVLST(SURF.idx_struct(:,i),3) = SURF.matNPVLST(SURF.idx_struct(:,i),3) + VEHI.matVEHUVW(:,3).*COND.valDELTIME + (OUTP.matDEFGLOB(valTIMESTEP,i) - OUTP.matDEFGLOB(valTIMESTEP-1,i)) + rotNPVLST(SURF.idx_struct(:,i),3);
    SURF.matNPVLST(SURF.idx_struct(:,i),1) = SURF.matNPVLST(SURF.idx_struct(:,i),1) + VEHI.matVEHUVW(:,1).*COND.valDELTIME + rotNPVLST(SURF.idx_struct(:,i),1);
    
    SURF.matEALST(SURF.idx_struct(:,i),3) = SURF.matEALST(SURF.idx_struct(:,i),3) + VEHI.matVEHUVW(:,3).*COND.valDELTIME + (OUTP.matDEFGLOB(valTIMESTEP,i) - OUTP.matDEFGLOB(valTIMESTEP-1,i));
    SURF.matEALST(SURF.idx_struct(:,i),1) = SURF.matEALST(SURF.idx_struct(:,i),1) + VEHI.matVEHUVW(:,1).*COND.valDELTIME;

    % Move DVE pts and elastic axis in y direction based on wing bending
    if i > 1
        SURF.matNPVLST(SURF.idx_struct(:,i),2) = SURF.matNPVLST(SURF.idx_struct(:,i-1),2) + (OUTP.matDEFGLOB(valTIMESTEP,i) - OUTP.matDEFGLOB(valTIMESTEP,i-1))/tan(OUTP.matSLOPE(valTIMESTEP,i));
        SURF.matEALST(SURF.idx_struct(:,i),2) = SURF.matEALST(SURF.idx_struct(:,i-1),2) + (OUTP.matDEFGLOB(valTIMESTEP,i) - OUTP.matDEFGLOB(valTIMESTEP,i-1))/tan(OUTP.matSLOPE(valTIMESTEP,i));        
    end
    
    % Average y-location of elastic axis to find y-location of wing CG
    if i > 1
        SURF.vecWINGCG(i-1,2) = (SURF.matEALST(SURF.idx_struct(1,i),2)+SURF.matEALST(SURF.idx_struct(1,i-1),2))/2;
    end
    
    translateNPVLST(SURF.idx_struct(:,i),1) = VEHI.matVEHUVW(:,1).*COND.valDELTIME;

end

if isempty(VEHI.vecPROPLOC) == 0
    VEHI.vecPROPLOC(:,3) = VEHI.vecPROPLOC(:,3) + VEHI.matVEHUVW(:,3).*COND.valDELTIME + (prop_def(valTIMESTEP) - prop_def(valTIMESTEP-1));
    VEHI.vecPROPLOC(:,1) = VEHI.vecPROPLOC(:,1) + VEHI.matVEHUVW(:,1).*COND.valDELTIME;
    
    deltaEPSPROP = (prop_twist(valTIMESTEP)-prop_twist(valTIMESTEP-1));
    ROTPROP = [cos(-deltaEPSPROP) 0 sin(-deltaEPSPROP); 0 1 0; -sin(-deltaEPSPROP) 0 cos(-deltaEPSPROP)];
    
    VEHI.vecPROPDIR = VEHI.vecPROPDIR*ROTPROP;
end

% Move tail pts and vehicle CG based on flight speed
SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matVEHUVW,size(SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:),1),1);

INPU.vecVEHCG = INPU.vecVEHCG + VEHI.matVEHUVW.*COND.valDELTIME;

% Rotate vehicle about CG based on flight-dynamics
if FLAG.FLIGHTDYN == 1
%     tempNPVLST = SURF.matNPVLST - INPU.vecVEHCG;
%     deltaEPS = (TRIM.perturb(end,4)-TRIM.perturb(end-1,4));
%     ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];
% 
%     vlst2 = (ROT*tempNPVLST')' + INPU.vecVEHCG;
% 
%     rotNPVLST2 = vlst2 - SURF.matNPVLST;
    [rotNPVLST2, ~] = fcnYROT(TRIM.perturb(end,4)-TRIM.perturb(end-1,4),SURF.matNPVLST,INPU.vecVEHCG);
    [~, VEHI.vecPAYLCG] = fcnYROT(TRIM.perturb(end,4)-TRIM.perturb(end-1,4),VEHI.vecPAYLCG,INPU.vecVEHCG);
    [~, VEHI.vecFUSECG] = fcnYROT(TRIM.perturb(end,4)-TRIM.perturb(end-1,4),VEHI.vecFUSECG,INPU.vecVEHCG);
    [~, VEHI.vecWINGCG(2,:)] = fcnYROT(TRIM.perturb(end,4)-TRIM.perturb(end-1,4),VEHI.vecWINGCG(2,:),INPU.vecVEHCG);
else
    rotNPVLST2 = zeros(size(SURF.matNPVLST,1),3);
end

VEHI.vecPAYLCG = VEHI.vecPAYLCG + VEHI.matVEHUVW.*COND.valDELTIME;
VEHI.vecFUSECG = VEHI.vecFUSECG + VEHI.matVEHUVW.*COND.valDELTIME;
VEHI.vecWINGCG(2,:) = VEHI.vecWINGCG(2,:) + VEHI.matVEHUVW.*COND.valDELTIME;
SURF.matNPVLST = SURF.matNPVLST + rotNPVLST2;

% update INPU.matVEHORIG positions
INPU.matVEHORIG = INPU.matVEHORIG + VEHI.matVEHUVW.*COND.valDELTIME;

% moveNPVLST = moveNPVLST + rotNPVLST2;

%% Move wing and generate new wake elements

% Old trailing edge vertices
MISC.matNEWWAKE(:,:,4) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,3) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),4),:);
MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),3),:);

% Update SURF.matVLST and SURF.matNTVLST
% SURF.matNTVLST = SURF.matNTVLST + rotNPVLST;

if any(FLAG.vecTRIMABLE == 1) == 1
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matVEHUVW*COND.valDELTIME;
    
end

% INPU.vecVEHCG = INPU.vecVEHCG + VEHI.matVEHUVW*COND.valDELTIME;

% if any(FLAG.vecFLEXIBLE == 1) == 1
%    
%     SURF.matEALST = SURF.matEALST + VEHI.matVEHUVW*COND.valDELTIME;
%     
% end

% Move stiff tail if it exists
% SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matVEHUVW,size(SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:),1),1);
% SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matVEHUVW,size(SURF.matNTVLST(tempENDWING+1:tempENDTAIL,:),1),1);

