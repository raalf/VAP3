function [SURF, INPU, MISC, VISC] = fcnMOVESURFACE_STATIC(INPU, VEHI, MISC, COND, SURF, VISC)


% Do the following need updating?
% VEHI.matVEHROT

INPU.valVEHICLES = length(VEHI.matVEHUVW(:,1));

% pre-calculate rad per timestep of rotors
vecROTORRADPS = COND.vecROTORRPM.*2.*pi./60;
vecROTORDEL = vecROTORRADPS.*COND.valDELTIME;
% TODO: insert warning if vecROTORRAD>2*pi

% update INPU.matVEHORIG positions
INPU.matVEHORIG = INPU.matVEHORIG + VEHI.matVEHUVW.*COND.valDELTIME;

% crate vecVLSTVEH which is a lookup vector for vertice to vehicle ID
vecVLSTVEH = unique([reshape(SURF.matDVE,[],1), repmat(SURF.vecDVEVEHICLE,4,1)],'rows');
vecVLSTVEH = vecVLSTVEH(:,2);

vecNTVLSTVEH = unique([reshape(SURF.matNPDVE,[],1), repmat(SURF.vecDVEVEHICLE,4,1)],'rows');
vecNTVLSTVEH = vecNTVLSTVEH(:,2);
% check error (should never fail if no vertices are shared between
% differnet vehicles
if length(vecVLSTVEH(:,1)) ~= length(SURF.matVLST(:,1))
    disp('vecVLSTVEH does not match vehicle ID');
end

% translation matrix for the vertice list
SURF.matVLSTTRANS = COND.valDELTIME.*VEHI.matVEHUVW(vecVLSTVEH,:);
SURF.matNTVLSTTRANS = COND.valDELTIME.*VEHI.matVEHUVW(vecNTVLSTVEH,:);
% translation matrix for the dve list
SURF.matDVETRANS  = COND.valDELTIME.*VEHI.matVEHUVW(SURF.vecDVEVEHICLE,:);
% [ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,~, ~, ~, ~, SURF.matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(SURF.matCENTER, SURF.matVLST(SURF.matDVE(:,1),:), SURF.matVLST(SURF.matDVE(:,2),:), SURF.matVLST(SURF.matDVE(:,3),:), SURF.matVLST(SURF.matDVE(:,4),:) );


% SURF.matDVETRANS holds UINF of each DVE due to tranlsation of vehicle
% hence excluding the effect of rotating rotors
SURF.matUINFVEH = -VEHI.matVEHUVW(SURF.vecDVEVEHICLE,:);
SURF.matUINFROT = SURF.matUINFVEH.*0;
SURF.matUINFTE = SURF.matUINFVEH.*0;

% Old trailing edge vertices
MISC.matNEWWAKE(:,:,4) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,3) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate WAKE.matWADJE)
MISC.matNPNEWWAKE(:,:,4) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,3) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);

% Translate Vehicles
SURF.matVLST = SURF.matVLST + SURF.matVLSTTRANS;
SURF.matCENTER = SURF.matCENTER + SURF.matDVETRANS;
SURF.matNTVLST = SURF.matNTVLST + SURF.matNTVLSTTRANS;
SURF.matNPVLST = SURF.matNPVLST + SURF.matNTVLSTTRANS;


% Circling Flight
% "backtrack" the UVW translation from previous lines of code, and apply the circling instead
for n = 1:INPU.valVEHICLES
   if ~isnan(VEHI.vecVEHRADIUS(n)) == 1
        idxDVEVEH = SURF.vecDVEVEHICLE == n;
        idxVLSTVEH = unique(SURF.matDVE(idxDVEVEH,:));
        idxFUSEVEH = VISC.vecFUSEVEHICLE == n;
        
        SURF.matVLST(idxVLSTVEH,1:2)   = SURF.matVLST(idxVLSTVEH,1:2)   - SURF.matVLSTTRANS(idxVLSTVEH,1:2);
        SURF.matNTVLST(idxVLSTVEH,1:2) = SURF.matNTVLST(idxVLSTVEH,1:2) - SURF.matVLSTTRANS(idxVLSTVEH,1:2);
        SURF.matCENTER(idxDVEVEH,1:2)  = SURF.matCENTER(idxDVEVEH,1:2)  - SURF.matDVETRANS(idxDVEVEH,1:2);
        VISC.matFUSEGEOM(:,:,1:2,idxFUSEVEH) = VISC.matFUSEGEOM(:,:,1:2,idxFUSEVEH) + matFUSETRANS(:,:,1:2,idxFUSEVEH);

        % reposition to vecCIRORIG to rotate the vehicle
        SURF.matVLST(idxVLSTVEH,1:2)   = SURF.matVLST(idxVLSTVEH,1:2)   - MISC.matCIRORIG(n,1:2);
        SURF.matNTVLST(idxVLSTVEH,1:2) = SURF.matNTVLST(idxVLSTVEH,1:2) - MISC.matCIRORIG(n,1:2);
        SURF.matCENTER(idxDVEVEH,1:2)  = SURF.matCENTER(idxDVEVEH,1:2)  - MISC.matCIRORIG(n,1:2);
%         VISC.matFUSEGEOM(:,:,1:2,idxFUSEVEH) = VISC.matFUSEGEOM(:,:,1:2,idxFUSEVEH) + MISC.matCIRORIG(n,1:2);
        
        % rotate(YAW) vehicle by VEHI.matVEHROTRATE(n,3)*COND.valDELTIME
        dcmVEHSTEP = angle2dcm(-VEHI.matVEHROTRATE(n,3)*COND.valDELTIME,0,0,'ZXY');
        
        SURF.matVLST(idxVLSTVEH,:)   = SURF.matVLST(idxVLSTVEH,:)   * dcmVEHSTEP;
        SURF.matNTVLST(idxVLSTVEH,:) = SURF.matNTVLST(idxVLSTVEH,:) * dcmVEHSTEP;
        SURF.matCENTER(idxDVEVEH,:)  = SURF.matCENTER(idxDVEVEH,:)  * dcmVEHSTEP;
      
        % 
        SURF.matVLST(idxVLSTVEH,1:2)   = SURF.matVLST(idxVLSTVEH,1:2)   + MISC.matCIRORIG(n,1:2);
        SURF.matNTVLST(idxVLSTVEH,1:2) = SURF.matNTVLST(idxVLSTVEH,1:2) + MISC.matCIRORIG(n,1:2);
        SURF.matCENTER(idxDVEVEH,1:2)  = SURF.matCENTER(idxDVEVEH,1:2)  + MISC.matCIRORIG(n,1:2);
%       VEHI.matVEHROTRATE(n,:)
%       MISC.matCIRORIG(n,:)
   end
end

% Rotate Rotors
valROTORS = length(VEHI.vecROTORVEH);
for n = 1:valROTORS
    
    dcmHUB2GLOB = angle2dcm(VEHI.matVEHROT(VEHI.vecROTORVEH(n),3),VEHI.matVEHROT(VEHI.vecROTORVEH(n),1),VEHI.matVEHROT(VEHI.vecROTORVEH(n),2),'ZXY');
%     dcmXY2HUB = quat2dcm(axang2quat(vrrotvec(INPU.VISC.matROTORAXIS(n,:),[0 0 1])));
    dcmXY2HUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], INPU.matROTORAXIS(n,:))));
    dcmROTORSTEP = angle2dcm(vecROTORDEL(n),0,0,'ZXY');

    % pre-calculate trans and rot matrices
    transGLOB2VEH = INPU.matROTORHUB(n,:) * dcmHUB2GLOB + INPU.matVEHORIG(VEHI.vecROTORVEH(n),:);
    
    
    idxDVEROTOR = SURF.vecDVEROTOR==n;
    idxVLSTROTOR = unique(SURF.matDVE(idxDVEROTOR,:));
    
    tempROTORVLST = SURF.matVLST(idxVLSTROTOR,:);
    tempROTORVLST = tempROTORVLST - transGLOB2VEH;
    
    tempROTORNTVLST = SURF.matNTVLST(idxVLSTROTOR,:);
    tempROTORNTVLST = tempROTORNTVLST - transGLOB2VEH;
           
    tempROTORCENTER = SURF.matCENTER(idxDVEROTOR,:);
    tempROTORCENTER = tempROTORCENTER - transGLOB2VEH;
           
        
    % transform rotor from global to hub plane
    tempROTORVLST = tempROTORVLST / dcmHUB2GLOB;
    tempROTORNTVLST = tempROTORNTVLST / dcmHUB2GLOB;
    tempROTORCENTER = tempROTORCENTER / dcmHUB2GLOB;    
    
    % transform rotor from hub plane to xy plane
    tempROTORVLST = tempROTORVLST / dcmXY2HUB;
    tempROTORNTVLST = tempROTORNTVLST / dcmXY2HUB;
    tempROTORCENTER = tempROTORCENTER / dcmXY2HUB;

    % timestep rotor in local XY hub plane
    tempROTORVLST = tempROTORVLST * dcmROTORSTEP;
    tempROTORNTVLST = tempROTORNTVLST * dcmROTORSTEP;
    tempROTORCENTER = tempROTORCENTER * dcmROTORSTEP;
    tempROTORUINF = cross(repmat([0,0,-vecROTORRADPS(n)],length(tempROTORCENTER(:,1)),1),tempROTORCENTER);    

    % transform rotor from xy plane to hub plane
    tempROTORVLST = tempROTORVLST * dcmXY2HUB;
    tempROTORNTVLST = tempROTORNTVLST * dcmXY2HUB;
    tempROTORCENTER = tempROTORCENTER * dcmXY2HUB;
    tempROTORUINF = tempROTORUINF * dcmXY2HUB;
    
    % transform rotor from hub plane to global
    tempROTORVLST = tempROTORVLST * dcmHUB2GLOB;
    tempROTORNTVLST = tempROTORNTVLST * dcmHUB2GLOB;
    tempROTORCENTER = tempROTORCENTER * dcmHUB2GLOB;
    tempROTORUINF = tempROTORUINF * dcmHUB2GLOB;
    
    % write rotated rotor to SURF.matVLST
    SURF.matVLST(idxVLSTROTOR,:) = tempROTORVLST + transGLOB2VEH;
    SURF.matNTVLST(idxVLSTROTOR,:) = tempROTORNTVLST + transGLOB2VEH;
    SURF.matCENTER(idxDVEROTOR,:) = tempROTORCENTER + transGLOB2VEH;
    SURF.matUINFROT(idxDVEROTOR,:) = tempROTORUINF;
end

% combine SURF.matUINFROTOR and SURF.matUINFVEH
SURF.matUINF = SURF.matUINFROT + SURF.matUINFVEH;

% [~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, ~, ~, ~, ~, SURF.matDVENORM, ~, ~, ~] = fcnVLST2DVEPARAM(SURF.matDVE, SURF.matVLST);
[ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,~, ~, ~, ~, SURF.matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(SURF.matCENTER, SURF.matVLST(SURF.matDVE(:,1),:), SURF.matVLST(SURF.matDVE(:,2),:), SURF.matVLST(SURF.matDVE(:,3),:), SURF.matVLST(SURF.matDVE(:,4),:) );

% New trailing edge vertices
MISC.matNEWWAKE(:,:,1) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,4),:);
MISC.matNEWWAKE(:,:,2) = SURF.matVLST(SURF.matDVE(SURF.vecDVETE>0,3),:);

% New non-planar trailing edge vertices (used to calculate WAKE.matWADJE)
MISC.matNPNEWWAKE(:,:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,4),:);
MISC.matNPNEWWAKE(:,:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE>0,3),:);
