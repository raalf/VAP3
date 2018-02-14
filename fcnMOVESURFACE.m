function [matUINF, matUINFTE, matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, ...
    matFUSEGEOM, matNEWWAKEPANEL, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matDVENORM, ...
    matNTVLST, matUINFROT] = fcnMOVESURFACE( matVEHORIG, matVEHUVW, ...
    matVEHROTRATE, matCIRORIG, vecVEHRADIUS, valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE, matFUSEGEOM, vecFUSEVEHICLE, ...
    matVEHROT, vecROTORVEH, matROTORHUB, matROTORAXIS, vecDVEROTOR, vecROTORRPM, matPANELTE, matNTVLST)


% Do the following need updating?
% matVEHROT


% Saving panel corner points, this helps with triangular wake generation
matNEWWAKEPANEL = zeros(size(matPANELTE,1),3,4);
matNEWWAKEPANEL(:,:,3) = matVLST(matPANELTE(:,2),:);
matNEWWAKEPANEL(:,:,4) = matVLST(matPANELTE(:,1),:);


valVEHICLES = length(matVEHUVW(:,1));

% pre-calculate rad per timestep of rotors
vecROTORRADPS = vecROTORRPM.*2.*pi./60;
vecROTORDEL = vecROTORRADPS.*valDELTIME;
% TODO: insert warning if vecROTORRAD>2*pi

% update matVEHORIG positions
matVEHORIG = matVEHORIG + matVEHUVW.*valDELTIME;

% crate vecVLSTVEH which is a lookup vector for vertice to vehicle ID
vecVLSTVEH = unique([reshape(matDVE,[],1), repmat(vecDVEVEHICLE,4,1)],'rows');
vecVLSTVEH = vecVLSTVEH(:,2);
% check error (should never fail if no vertices are shared between
% differnet vehicles
if length(vecVLSTVEH(:,1)) ~= length(matVLST(:,1))
    disp('vecVLSTVEH does not match vehicle ID');
end

% translation matrix for the vertice list
matVLSTTRANS = valDELTIME.*matVEHUVW(vecVLSTVEH,:);
% translation matrix for the dve list
matDVETRANS  = valDELTIME.*matVEHUVW(vecDVEVEHICLE,:);
% [ ~, ~, vecDVEROLL, vecDVEPITCH, vecDVEYAW,~, ~, ~, ~, matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(matCENTER, matVLST(matDVE(:,1),:), matVLST(matDVE(:,2),:), matVLST(matDVE(:,3),:), matVLST(matDVE(:,4),:) );



% Fuselage
matFUSETRANS = valDELTIME.*matVEHUVW(vecFUSEVEHICLE,:);
sz = size(matFUSEGEOM);
matFUSETRANS = repmat(reshape(matFUSETRANS',1,1,3,length(vecFUSEVEHICLE)),sz(1),sz(2),1,1);
matFUSEGEOM = matFUSEGEOM + matFUSETRANS;



% matDVETRANS holds UINF of each DVE due to tranlsation of vehicle
% hence excluding the effect of rotating rotors
matUINFVEH = -matVEHUVW(vecDVEVEHICLE,:);
matUINFROT = matUINFVEH.*0;
matUINFTE = matUINFVEH.*0;

% Old trailing edge vertices
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,4) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNTVLST(matDVE(vecDVETE>0,3),:);

% Translate Vehicles
matVLST = matVLST + matVLSTTRANS;
matCENTER = matCENTER + matDVETRANS;
matNTVLST = matNTVLST + matVLSTTRANS;


% Circling Flight
% "backtrack" the UVW translation from previous lines of code, and apply the circling instead
for n = 1:valVEHICLES
   if ~isnan(vecVEHRADIUS(n)) == 1
        idxDVEVEH = vecDVEVEHICLE == n;
        idxVLSTVEH = unique(matDVE(idxDVEVEH,:));
        idxFUSEVEH = vecFUSEVEHICLE == n;
        
        matVLST(idxVLSTVEH,1:2)   = matVLST(idxVLSTVEH,1:2)   - matVLSTTRANS(idxVLSTVEH,1:2);
        matNTVLST(idxVLSTVEH,1:2) = matNTVLST(idxVLSTVEH,1:2) - matVLSTTRANS(idxVLSTVEH,1:2);
        matCENTER(idxDVEVEH,1:2)  = matCENTER(idxDVEVEH,1:2)  - matDVETRANS(idxDVEVEH,1:2);
        matFUSEGEOM(:,:,1:2,idxFUSEVEH) = matFUSEGEOM(:,:,1:2,idxFUSEVEH) + matFUSETRANS(:,:,1:2,idxFUSEVEH);

        % reposition to vecCIRORIG to rotate the vehicle
        matVLST(idxVLSTVEH,1:2)   = matVLST(idxVLSTVEH,1:2)   - matCIRORIG(n,1:2);
        matNTVLST(idxVLSTVEH,1:2) = matNTVLST(idxVLSTVEH,1:2) - matCIRORIG(n,1:2);
        matCENTER(idxDVEVEH,1:2)  = matCENTER(idxDVEVEH,1:2)  - matCIRORIG(n,1:2);
%         matFUSEGEOM(:,:,1:2,idxFUSEVEH) = matFUSEGEOM(:,:,1:2,idxFUSEVEH) + matCIRORIG(n,1:2);
        
        % rotate(YAW) vehicle by matVEHROTRATE(n,3)*valDELTIME
        dcmVEHSTEP = angle2dcm(-matVEHROTRATE(n,3)*valDELTIME,0,0,'ZXY');
        
        matVLST(idxVLSTVEH,:)   = matVLST(idxVLSTVEH,:)   * dcmVEHSTEP;
        matNTVLST(idxVLSTVEH,:) = matNTVLST(idxVLSTVEH,:) * dcmVEHSTEP;
        matCENTER(idxDVEVEH,:)  = matCENTER(idxDVEVEH,:)  * dcmVEHSTEP;
      
        % 
        matVLST(idxVLSTVEH,1:2)   = matVLST(idxVLSTVEH,1:2)   + matCIRORIG(n,1:2);
        matNTVLST(idxVLSTVEH,1:2) = matNTVLST(idxVLSTVEH,1:2) + matCIRORIG(n,1:2);
        matCENTER(idxDVEVEH,1:2)  = matCENTER(idxDVEVEH,1:2)  + matCIRORIG(n,1:2);
%       matVEHROTRATE(n,:)
%       matCIRORIG(n,:)
   end
end

% Rotate Rotors
valROTORS = length(vecROTORVEH);
for n = 1:valROTORS
    
    dcmHUB2GLOB = angle2dcm(matVEHROT(vecROTORVEH(n),3),matVEHROT(vecROTORVEH(n),1),matVEHROT(vecROTORVEH(n),2),'ZXY');
%     dcmXY2HUB = quat2dcm(axang2quat(vrrotvec(matROTORAXIS(n,:),[0 0 1])));
    dcmXY2HUB = quat2dcm(axang2quat(vrrotvec([0 0 1], matROTORAXIS(n,:))));
    dcmROTORSTEP = angle2dcm(vecROTORDEL(n),0,0,'ZXY');

    % pre-calculate trans and rot matrices
    transGLOB2VEH = matROTORHUB(n,:) * dcmHUB2GLOB + matVEHORIG(vecROTORVEH(n),:);
    
    
    idxDVEROTOR = vecDVEROTOR==n;
    idxVLSTROTOR = unique(matDVE(idxDVEROTOR,:));
    
    tempROTORVLST = matVLST(idxVLSTROTOR,:);
    tempROTORVLST = tempROTORVLST - transGLOB2VEH;
    
    tempROTORNTVLST = matNTVLST(idxVLSTROTOR,:);
    tempROTORNTVLST = tempROTORNTVLST - transGLOB2VEH;
           
    tempROTORCENTER = matCENTER(idxDVEROTOR,:);
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
    
    % write rotated rotor to matVLST
    matVLST(idxVLSTROTOR,:) = tempROTORVLST + transGLOB2VEH;
    matNTVLST(idxVLSTROTOR,:) = tempROTORNTVLST + transGLOB2VEH;
    matCENTER(idxDVEROTOR,:) = tempROTORCENTER + transGLOB2VEH;
    matUINFROT(idxDVEROTOR,:) = tempROTORUINF;
end

% combine matUINFROTOR and matUINFVEH
matUINF = matUINFROT + matUINFVEH;

% [~, ~, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ~, ~, ~, ~, matDVENORM, ~, ~, ~] = fcnVLST2DVEPARAM(matDVE, matVLST);
[ ~, ~, vecDVEROLL, vecDVEPITCH, vecDVEYAW,~, ~, ~, ~, matDVENORM, ~, ~, ~, ~] = fcnDVECORNER2PARAM(matCENTER, matVLST(matDVE(:,1),:), matVLST(matDVE(:,2),:), matVLST(matDVE(:,3),:), matVLST(matDVE(:,4),:) );

% New trailing edge vertices
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNTVLST(matDVE(vecDVETE>0,3),:);

matNEWWAKEPANEL(:,:,1) = matVLST(matPANELTE(:,1),:);
matNEWWAKEPANEL(:,:,2) = matVLST(matPANELTE(:,2),:);
