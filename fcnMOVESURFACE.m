function [matUINF, matUINFTE, matVEHORIG, ...
    matVLST, matCENTER, ...
    matNEWWAKE, matNPNEWWAKE, ...
    matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, ...
    valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, ...
    vecDVETE, matNTVLST, matFUSEGEOM, vecFUSEVEHICLE, ...
    matVEHROT, vecROTORVEH, matROTORHUBGLOB, ...
    matROTORHUB, matROTORAXIS, vecDVEROTOR, vecROTORRPM)
% This function moves a wing (NOT rotor) by translating all of the vertices
% in the VLST and the in-centers of each triangle in CENTER.

% INPUT:
%   matVEHORIG - Current vehicle origin of each vehicle valVEHICLE x 3
%   matVEHUVW - velocity components of each vehicle valVEHICLE x 4
%   valDELTIME - Timestep size for this move
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
%   matELST - List of edges from fcnTRIANG
%   vecTE - List of trailing edge edges
% OUTPUT:
%   matUINF - Local 
%   matVEHORIG - vehicle origin after timestep of each vehicle valVEHICLE x 3
%   matVLST - New vertex list with moved points
%   matCENTER - New in-center list with moved points
%   matNEWWAKE - Outputs a 4 x 3 x n matrix of points for the wake DVE generation


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



% Fuselage
matFUSETRANS = valDELTIME.*matVEHUVW(vecFUSEVEHICLE,:);
sz = size(matFUSEGEOM);
matFUSETRANS = repmat(reshape(matFUSETRANS',1,1,3,length(vecFUSEVEHICLE)),sz(1),sz(2),1,1);
matFUSEGEOM = matFUSEGEOM + matFUSETRANS;








% matDVETRANS holds UINF of each DVE due to tranlsation of vehicle
% hence excluding the effect of rotating rotors
matUINFVEH = -matVEHUVW(vecDVEVEHICLE,:);
matUINFROTOR = matUINFVEH.*0;
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

% Rotate Rotors
valROTORS = length(vecROTORVEH);
for n = 1:valROTORS
    
    % pre-calculate trans and rot matrices
    transGLOB2VEH = matROTORHUBGLOB(n,:) + matVEHORIG(vecROTORVEH(n),:);
    dcmHUB2GLOB = angle2dcm(matVEHROT(vecROTORVEH(n),1),matVEHROT(vecROTORVEH(n),2),matVEHROT(vecROTORVEH(n),3),'XYZ');
    dcmXY2HUB = quat2dcm(axang2quat(vrrotvec(matROTORAXIS(n,:),[0 0 1])));
    dcmROTORSTEP = angle2dcm(0,0,vecROTORDEL(n),'XYZ');
    
    idxDVEROTOR = vecDVEROTOR==n;
    idxVLSTROTOR = unique(matDVE(idxDVEROTOR,:));
    
    tempROTORVLST = matVLST(idxVLSTROTOR,:);
    tempROTORVLST = tempROTORVLST - transGLOB2VEH;
       
    tempROTORCENTER = matCENTER(idxDVEROTOR,:);
    tempROTORCENTER = tempROTORCENTER - transGLOB2VEH;
           
        
    % transform rotor from global to hub plane
    tempROTORVLST = tempROTORVLST / dcmHUB2GLOB;
    tempROTORCENTER = tempROTORCENTER / dcmHUB2GLOB;    
    
    % transform rotor from hub plane to xy plane
    tempROTORVLST = tempROTORVLST / dcmXY2HUB;
    tempROTORCENTER = tempROTORCENTER / dcmXY2HUB;

    % timestep rotor in local XY hub plane
    tempROTORVLST = tempROTORVLST * dcmROTORSTEP;
    tempROTORCENTER = tempROTORCENTER * dcmROTORSTEP;
    tempROTORUINF = cross(repmat([0,0,-vecROTORRADPS(n)],length(tempROTORCENTER(:,1)),1),tempROTORCENTER);    
    
    % transform rotor from xy plane to hub plane
    tempROTORVLST = tempROTORVLST * dcmXY2HUB;
    tempROTORCENTER = tempROTORCENTER * dcmXY2HUB;
    tempROTORUINF = tempROTORUINF * dcmXY2HUB;
    
    % transform rotor from hub plane to global
    tempROTORVLST = tempROTORVLST * dcmHUB2GLOB;
    tempROTORCENTER = tempROTORCENTER * dcmHUB2GLOB;
    tempROTORUINF = tempROTORUINF * dcmHUB2GLOB;
    
    % write rotated rotor to matVLST
    matVLST(idxVLSTROTOR,:) = tempROTORVLST + transGLOB2VEH;
    matCENTER(idxDVEROTOR,:) = tempROTORCENTER + transGLOB2VEH;
    matUINFROTOR(idxDVEROTOR,:) = tempROTORUINF;
end


% combine matUINFROTOR and matUINFVEH
matUINF = matUINFROTOR + matUINFVEH;

% New trailing edge vertices
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNTVLST(matDVE(vecDVETE>0,3),:);