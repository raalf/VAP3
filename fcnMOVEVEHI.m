function [SURF, MISC, COND, INPU] = fcnMOVEVEHI(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM)

%% Move a vehicle based on the flight dynamic response

% Rotate vehicle about CG based on pitch response
% tempNPVLST = SURF.matNPVLST - INPU.vecVEHCG;
% 
% if FLAG.FLIGHTDYN == 1
%     deltaEPS = (TRIM.perturb(end,4)-TRIM.perturb(end-1,4));
%     ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];
% 
%     vlst2 = (ROT*tempNPVLST')' + INPU.vecVEHCG;
% 
%     rotNPVLST = vlst2 - SURF.matNPVLST;
% else
%     rotNPVLST = [];
% end

% glob to local translation
matNPVLST = SURF.matNPVLST - repmat(INPU.vecVEHCG,size(SURF.matNPVLST,1),1);
matVLST = SURF.matVLST - repmat(INPU.vecVEHCG,size(SURF.matVLST,1),1);
matCENTER = SURF.matCENTER - repmat(INPU.vecVEHCG,size(SURF.matCENTER,1),1);

% rotate
%     dcm = angle2dcm(matVEHROT(n,1), matVEHROT(n,2), matVEHROT(n,3), 'XYZ');
matVEHROT = [0, (TRIM.perturb(end,4)-TRIM.perturb(end-1,4)), 0];
dcm = angle2dcm(matVEHROT(3), matVEHROT(1), matVEHROT(2), 'ZXY');
matNPVLST = matNPVLST*dcm;
matVLST = matVLST*dcm;
matCENTER = matCENTER*dcm;

% local to global translation
matNPVLST = matNPVLST + repmat(INPU.vecVEHCG,size(SURF.matNPVLST,1),1);
matVLST = matVLST + repmat(INPU.vecVEHCG,size(SURF.matVLST,1),1);
matCENTER = matCENTER + repmat(INPU.vecVEHCG,size(SURF.matCENTER,1),1);

rotNPVLST = matNPVLST - SURF.matNPVLST;

% Translate vehicle based on freestream velocity and translational dynamics
% translateNPVLST(:,3) = VEHI.matVEHUVW(:,3).*COND.valDELTIME;
% translateNPVLST(:,2) = VEHI.matVEHUVW(:,2).*COND.valDELTIME;
% translateNPVLST(:,1) = VEHI.matVEHUVW(:,1).*COND.valDELTIME;

% moveNPVLST = translateNPVLST + rotNPVLST;

%% Move wing and generate new wake elements

% Old trailing edge vertices
MISC.matNEWWAKE(:,:,4) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,3) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(:,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);

% Update SURF.matVLST and SURF.matNTVLST
% SURF.matNPVLST = SURF.matNPVLST + moveNPVLST;
SURF.matNPVLST = matNPVLST + repmat(COND.valDELTIME.*VEHI.matVEHUVW,size(matNPVLST,1),1);
SURF.matVLST = matVLST + repmat(COND.valDELTIME.*VEHI.matVEHUVW,size(matVLST,1),1);
SURF.matCENTER = matCENTER + repmat(COND.valDELTIME.*VEHI.matVEHUVW,size(matCENTER,1),1);

INPU.matVEHORIG = INPU.matVEHORIG + VEHI.matVEHUVW.*COND.valDELTIME;
INPU.vecVEHCG = INPU.vecVEHCG + VEHI.matVEHUVW.*COND.valDELTIME;

[ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,~, ~, ~, ~, SURF.matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(SURF.matCENTER, SURF.matVLST(SURF.matDVE(:,1),:), SURF.matVLST(SURF.matDVE(:,2),:), SURF.matVLST(SURF.matDVE(:,3),:), SURF.matVLST(SURF.matDVE(:,4),:) );

if any(FLAG.vecTRIMABLE == 1) == 1
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matVEHUVW*COND.valDELTIME;
    
end

% [ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
%     SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
%     = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);

% New non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(:,:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);

% New trailing edge vertices
MISC.matNEWWAKE(:,:,1) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,2) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);
