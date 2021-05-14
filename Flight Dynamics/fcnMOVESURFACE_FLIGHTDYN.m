function [SURF, MISC, COND, INPU, VEHI] = fcnMOVESURFACE_FLIGHTDYN(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP)

tempENDWING = max(max(SURF.matNPDVE(SURF.idxFLEX,:))); % Index for last row in matNTVLST that corresponds to the wing
tempENDTAIL = max(max(SURF.matNPDVE(SURF.idxTAIL,:)));
    
% Translation matrix for NPVLST points
for i = 1:size(SURF.idx_struct,2)
    
    SURF.matEALST(SURF.idx_struct(:,i),3) = SURF.matEALST(SURF.idx_struct(:,i),3) + VEHI.matGLOBUVW(:,3).*COND.valDELTIME;
    SURF.matEALST(SURF.idx_struct(:,i),1) = SURF.matEALST(SURF.idx_struct(:,i),1) + VEHI.matGLOBUVW(:,1).*COND.valDELTIME;
    
    % Average y-location of elastic axis to find y-location of wing CG
    if i > 1
        SURF.vecWINGCG(i-1,2) = (SURF.matEALST(SURF.idx_struct(1,i),2)+SURF.matEALST(SURF.idx_struct(1,i-1),2))/2;
    end

end

% Rotate vehicle about CG based on flight-dynamics
if FLAG.FLIGHTDYN == 1
    [rotNPVLST2, ~] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.matNPVLST,INPU.vecVEHCG);
    [~, VEHI.vecPAYLCG] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecPAYLCG,INPU.vecVEHCG);
    [~, VEHI.vecFUSECG] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecFUSECG,INPU.vecVEHCG);
%     [~, VEHI.vecWINGCG(2,:)] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecWINGCG(2,:),INPU.vecVEHCG);
    [~, VEHI.vecWINGCG] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecWINGCG,INPU.vecVEHCG);
else
    rotNPVLST2 = zeros(size(SURF.matNPVLST,1),3);
end

INPU.vecVEHCG = INPU.vecVEHCG + VEHI.matGLOBUVW.*COND.valDELTIME;
SURF.matEA = SURF.matEA + VEHI.matGLOBUVW*COND.valDELTIME;
SURF.matCG = SURF.matCG + VEHI.matGLOBUVW*COND.valDELTIME;
VEHI.vecPAYLCG = VEHI.vecPAYLCG + VEHI.matGLOBUVW.*COND.valDELTIME;
VEHI.vecFUSECG = VEHI.vecFUSECG + VEHI.matGLOBUVW.*COND.valDELTIME;
SURF.vecWINGCG = SURF.vecWINGCG + VEHI.matGLOBUVW.*COND.valDELTIME;
% VEHI.vecWINGCG(2,:) = VEHI.vecWINGCG(2,:) + VEHI.matGLOBUVW.*COND.valDELTIME;
VEHI.vecWINGCG = VEHI.vecWINGCG + VEHI.matGLOBUVW.*COND.valDELTIME;

SURF.matEALST(1:tempENDWING,:) = SURF.matEALST(1:tempENDWING,:) + VEHI.matGLOBUVW.*COND.valDELTIME;
SURF.matCGLST(1:tempENDWING,:) = SURF.matCGLST(1:tempENDWING,:) + VEHI.matGLOBUVW.*COND.valDELTIME;
SURF.matBEAMLOC(:,:,valTIMESTEP) = SURF.matEALST(1:INPU.valNSELE,:);
SURF.matBEAMCGLOC(:,:,valTIMESTEP) = SURF.matCGLST(1:INPU.valNSELE,:);
SURF.matAEROCNTR = SURF.matAEROCNTR + VEHI.matGLOBUVW.*COND.valDELTIME;

% update INPU.matVEHORIG positions
INPU.matVEHORIG = INPU.matVEHORIG + VEHI.matGLOBUVW.*COND.valDELTIME;

%% Move wing and generate new wake elements

% Old trailing edge vertices
MISC.matNEWWAKE(:,:,4) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,3) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(:,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);

SURF.matNPVLST = SURF.matNPVLST + rotNPVLST2;

% Move DVE pts and elastic axis in x and z direction based on wing
% bending and twist
SURF.matNPVLST(:,3) = SURF.matNPVLST(:,3) + VEHI.matGLOBUVW(:,3).*COND.valDELTIME;
SURF.matNPVLST(:,1) = SURF.matNPVLST(:,1) + VEHI.matGLOBUVW(:,1).*COND.valDELTIME;

MISC.matNPNEWWAKE(:,:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);
% MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),4),:);
% MISC.matNPNEWWAKE(length(find(SURF.vecDVETE(SURF.idxFLEX) == 3))+1:end,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.idxTAIL(SURF.vecDVETE(SURF.idxTAIL)>0),3),:);

% Update SURF.matVLST and SURF.matNTVLST
% SURF.matNTVLST = SURF.matNTVLST + rotNPVLST;

if any(FLAG.vecTRIMABLE == 1) == 1
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matGLOBUVW*COND.valDELTIME;
    
end