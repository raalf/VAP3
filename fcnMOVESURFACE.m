function [matUINF, matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE, matNTVLST, matFUSEGEOM, vecFUSEVEHICLE)
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

matFUSETRANS = valDELTIME.*matVEHUVW(vecFUSEVEHICLE,:);
sz = size(matFUSEGEOM);
matFUSETRANS = repmat(reshape(matFUSETRANS',1,1,3,length(vecFUSEVEHICLE)),sz(1),sz(2),1,1);
matFUSEGEOM = matFUSEGEOM + matFUSETRANS;

% matDVETRANS holds UINF of each DVE due to tranlsation of vehicle
% hence excluding the effect of rotating rotors
matUINF = matDVETRANS;


% Old trailing edge vertices
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,4) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNTVLST(matDVE(vecDVETE>0,3),:);

matVLST = matVLST + matVLSTTRANS;
matCENTER = matCENTER + matDVETRANS;
matNTVLST = matNTVLST + matVLSTTRANS;

% New trailing edge vertices
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNTVLST(matDVE(vecDVETE>0,3),:);