function [SURF, MISC, COND, INPU, VEHI, OUTP] = fcnMOVEFLEXVEHI2(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP)

tempENDWING = max(max(SURF.matNPDVE(SURF.idxFLEX,:))); % Index for last row in matNPVLST that corresponds to the wing
tempENDTAIL = max(max(SURF.matNPDVE(SURF.idxTAIL,:))); % Index for last row in matNPVLST that corresponds to the tail

% Deflection and twist at the wing control point y-locations
cpt_def(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP,:),SURF.center_dist);
cpt_twist(valTIMESTEP,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP,:),SURF.center_dist);

cpt_def(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matDEFGLOB(valTIMESTEP-1,:),SURF.center_dist);
cpt_twist(valTIMESTEP-1,:) = interp1(SURF.vecSPANLOC,OUTP.matTWISTGLOB(valTIMESTEP-1,:),SURF.center_dist);

tempEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),SURF.vecWINGCG(:,2));
tempEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.vecWINGCG(:,2));
tempEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),SURF.vecWINGCG(:,2));

% Translation matrix for NPVLST points
for i = 1:size(SURF.idx_struct,2)

    % Perform rotation due to wing twist about elastic axis for DVE pts
    tempNPVLST = SURF.matNPVLST(SURF.idx_struct(:,i),:) - SURF.matEALST(SURF.idx_struct(:,i),:);

    deltaEPS = (OUTP.matTWISTGLOB(valTIMESTEP,i)-OUTP.matTWISTGLOB(valTIMESTEP-1,i));
    ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];
    
    vlst2 = (ROT*tempNPVLST')' + SURF.matEALST(SURF.idx_struct(:,i),:);
    
    rotNPVLST(SURF.idx_struct(:,i),:) = vlst2 - SURF.matNPVLST(SURF.idx_struct(:,i),:);
    
    % Perform rotation due to wing twist about elastic axis for wing cg
    % pts
    if i > 1
        tempCG = (SURF.vecWINGCG(i-1,:) - tempEA(i-1,:))*0;

        deltaEPSCG = (cpt_twist(valTIMESTEP,i-1)-cpt_twist(valTIMESTEP-1,i-1));
        ROTCG = [cos(deltaEPSCG) 0 sin(deltaEPSCG); 0 1 0; -sin(deltaEPSCG) 0 cos(deltaEPSCG)];

        cg2 = (ROTCG*tempCG')' + SURF.vecWINGCG(i-1,:);

        rotCG(i-1,:) = cg2 - SURF.vecWINGCG(i-1,:);
        
        SURF.vecWINGCG(i-1,3) = SURF.vecWINGCG(i-1,3) + VEHI.matGLOBUVW(:,3).*COND.valDELTIME + (cpt_def(valTIMESTEP,i-1) - cpt_def(valTIMESTEP-1,i-1)) + rotCG(i-1,3);
        SURF.vecWINGCG(i-1,1) = SURF.vecWINGCG(i-1,1) + VEHI.matGLOBUVW(:,1).*COND.valDELTIME + rotCG(i-1,1);
    end
    
end
    
    rij = 2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecDVEWING == 1) == 1);

    for iii = 1:size(OUTP.matDEFGLOB,2)

        if iii == 1
            ry(1,iii) = 0;
            rz(1,iii) = OUTP.matDEFGLOB(valTIMESTEP,iii);
        else
            ry(1,iii) = rij(iii-1)*cos(asin((OUTP.matDEFGLOB(valTIMESTEP,iii)-OUTP.matDEFGLOB(valTIMESTEP,iii-1))/rij(iii-1)));
            ry2(1,iii) = rij(iii-1)*cos(asin((OUTP.matDEFGLOB(valTIMESTEP-1,iii)-OUTP.matDEFGLOB(valTIMESTEP-1,iii-1))/rij(iii-1)));
            rz(1,iii) = OUTP.matDEFGLOB(valTIMESTEP,iii)-OUTP.matDEFGLOB(valTIMESTEP-1,iii);
        end

    end
    
    for i = 1:size(SURF.idx_struct,2)
        SURF.matBEAMROLL(:,valTIMESTEP) = SURF.vecDVEROLL(1:INPU.valNSELE-1);
        dRy(unique(SURF.idx_struct(:,i)),1) = repmat(cumsum(ry(i))',size(unique(SURF.idx_struct(:,i)),1),1) - SURF.matNPVLST(unique(SURF.idx_struct(:,i)),2);
        dRy2(unique(SURF.idx_struct(:,i)),1) = repmat(cumsum(ry2(i))',size(unique(SURF.idx_struct(:,i)),1),1) - SURF.matNPVLST(unique(SURF.idx_struct(:,i)),2);
        dRy(unique(SURF.idx_struct(:,i)),1) = dRy(unique(SURF.idx_struct(:,i)),1) - dRy2(unique(SURF.idx_struct(:,i)),1);
        dRz(unique(SURF.idx_struct(:,i)),1) = repmat(rz(i)',size(unique(SURF.idx_struct(:,i)),1),1);
    end
    
    % Span vector of the beam
%     s = repmat([0 -1 0; -(SURF.matBEAMLOC(2:end,:,valTIMESTEP-1) - SURF.matBEAMLOC(1:end-1,:,valTIMESTEP-1))./rij],INPU.vecM(1)+1,1);
              
    elastic_translation = rotNPVLST + [zeros(size(dRz,1),1),dRy,zeros(size(dRz,1),1)] + [zeros(size(dRz,1),2),dRz];
%     elastic_translation = rotNPVLST + dRy.*s + [zeros(size(dRz,1),2),dRz];
    elastic_translation = fcnSTARGLOB(elastic_translation, deg2rad(COND.vecVEHROLL*ones(size(elastic_translation,1),1)), deg2rad(COND.vecVEHPITCH*ones(size(elastic_translation,1),1)), deg2rad(COND.vecVEHBETA*ones(size(elastic_translation,1),1)));
    
    % Move DVE pts and elastic axis based on wing bending and twist   
    SURF.matNPVLST(1:tempENDWING,:) = SURF.matNPVLST(1:tempENDWING,:) + VEHI.matGLOBUVW.*COND.valDELTIME + elastic_translation;
    
    SURF.matEALST(1:tempENDWING,:) = SURF.matEALST(1:tempENDWING,:) + VEHI.matGLOBUVW.*COND.valDELTIME + elastic_translation;
    SURF.matCGLST(1:tempENDWING,:) = SURF.matCGLST(1:tempENDWING,:) + VEHI.matGLOBUVW.*COND.valDELTIME + elastic_translation;
    SURF.matAEROCNTR = SURF.matAEROCNTR + VEHI.matGLOBUVW.*COND.valDELTIME + elastic_translation(1:sum(INPU.vecN(INPU.vecPANELWING == 1),1)+1,:);
    
    SURF.matBEAMOMEGA(:,:,valTIMESTEP) = [(3*SURF.matBEAMROLL(:,valTIMESTEP) - 4*SURF.matBEAMROLL(:,valTIMESTEP-1) + SURF.matBEAMROLL(:,valTIMESTEP-2))./(2*COND.valDELTIME),zeros(INPU.valNSELE-1,2)];    
    
    % Average y-location of elastic axis to find y-location of wing CG
    SURF.vecWINGCG(:,2) = (SURF.matEALST(2:SURF.idx_struct(1,end),2)+SURF.matEALST(1:SURF.idx_struct(1,end)-1,2))./2;

    % Move tail pts and vehicle CG based on vehicle velocity in global frame
    SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) = SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:) + COND.valDELTIME.*repmat(VEHI.matGLOBUVW,size(SURF.matNPVLST(tempENDWING+1:tempENDTAIL,:),1),1);

% Rotate vehicle about CG based on flight-dynamics
if FLAG.FLIGHTDYN == 1
    [rotNPVLST2, ~] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.matNPVLST,INPU.vecVEHCG);
    [~, VEHI.vecPAYLCG] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecPAYLCG,INPU.vecVEHCG);
    [~, VEHI.vecFUSEMASSLOC] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecFUSEMASSLOC,INPU.vecVEHCG);
    [~, SURF.vecWINGCG] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.vecWINGCG,INPU.vecVEHCG);
    [~, VEHI.vecWINGCG(2,:)] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),VEHI.vecWINGCG(2,:),INPU.vecVEHCG);
    [~, SURF.matAEROCNTR] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.matAEROCNTR,INPU.vecVEHCG);
    [~, SURF.matEALST(1:tempENDWING,:)] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.matEALST(1:tempENDWING,:),INPU.vecVEHCG);
    [~, SURF.matCGLST(1:tempENDWING,:)] = fcnYROT(VEHI.vecVEHDYN(valTIMESTEP,4)-VEHI.vecVEHDYN(valTIMESTEP-1,4),SURF.matCGLST(1:tempENDWING,:),INPU.vecVEHCG);
else
    rotNPVLST2 = zeros(size(SURF.matNPVLST,1),3);
end

SURF.matBEAMLOC(:,:,valTIMESTEP) = SURF.matEALST(1:INPU.valNSELE,:);
SURF.matBEAMCGLOC(:,:,valTIMESTEP) = SURF.matCGLST(1:INPU.valNSELE,:);

SURF.matBEAMVEL(:,:,valTIMESTEP) = (SURF.matBEAMLOC(:,:,valTIMESTEP) - SURF.matBEAMLOC(:,:,valTIMESTEP-1))./(COND.valDELTIME);

% Update position of all constants using vehicle velocity in global frame
INPU.vecVEHCG = INPU.vecVEHCG + VEHI.matGLOBUVW.*COND.valDELTIME;
VEHI.vecPAYLCG = VEHI.vecPAYLCG + VEHI.matGLOBUVW.*COND.valDELTIME;
VEHI.vecFUSEMASSLOC = VEHI.vecFUSEMASSLOC + VEHI.matGLOBUVW.*COND.valDELTIME;
VEHI.vecWINGCG(2,:) = VEHI.vecWINGCG(2,:) + VEHI.matGLOBUVW.*COND.valDELTIME;
SURF.matNPVLST = SURF.matNPVLST + rotNPVLST2;
INPU.matVEHORIG = INPU.matVEHORIG + VEHI.matGLOBUVW.*COND.valDELTIME;

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
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matGLOBUVW*COND.valDELTIME;
    
end

