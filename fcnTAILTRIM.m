function [SURF] = fcnTAILTRIM(SURF, FLAG, COND, deltaEPS, i)
% Function to twist the tail based on computed incidence angle

idxtail = SURF.vecDVEWING == 2; % DVEs that are on the tail

tempVLST = SURF.matVLST(unique(SURF.matDVE(idxtail,:)),:) - repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matVLST(unique(SURF.matDVE(idxtail,:)),1)),1);
tempNPVLST = SURF.matNPVLST(unique(SURF.matNPDVE(idxtail,:)),:) - repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matNPVLST(unique(SURF.matNPDVE(idxtail,:)),1)),1);
tempCENTER = SURF.matCENTER(idxtail,:) - repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matCENTER(idxtail,1)),1);

ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];

vlst2 = (ROT*tempVLST')' + repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matVLST(unique(SURF.matDVE(idxtail,:)),1)),1) ;
npvlst2 = (ROT*tempNPVLST')' + repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matNPVLST(unique(SURF.matNPDVE(idxtail,:)),1)),1) ;
center2 = (ROT*tempCENTER')' + repmat(SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:), length(SURF.matCENTER(idxtail,1)),1);

SURF.matVLST(unique(SURF.matDVE(idxtail,:)),:) = vlst2;
SURF.matNPVLST(unique(SURF.matNPDVE(idxtail,:)),:) = npvlst2;
SURF.matCENTER(idxtail,:) = center2;

% Store old tail pitch
SURF.vecTAILPITCHold = SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180;

[ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,~, ~, ~, ~, SURF.matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(SURF.matCENTER,...
    SURF.matVLST(SURF.matDVE(:,1),:), SURF.matVLST(SURF.matDVE(:,2),:), SURF.matVLST(SURF.matDVE(:,3),:), SURF.matVLST(SURF.matDVE(:,4),:) );

end