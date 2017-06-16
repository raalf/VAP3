function [matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE, matNTVLST)
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
%   matVEHORIG - vehicle origin after timestep of each vehicle valVEHICLE x 3
%   matVLST - New vertex list with moved points
%   matCENTER - New in-center list with moved points
%   matNEWWAKE - Outputs a 4 x 3 x n matrix of points for the wake DVE generation


matVEHORIG = matVEHORIG + matVEHUVW.*valDELTIME;
vecVLSTVEH = unique([reshape(matDVE,[],1), repmat(vecDVEVEHICLE,4,1)],'rows');
vecVLSTVEH = vecVLSTVEH(:,2);
if length(vecVLSTVEH(:,1)) ~= length(matVLST(:,1))
    disp('vecVLSTVEH does not match vehicle ID');
end



matVLSTTRANS = valDELTIME.*matVEHUVW(vecVLSTVEH,:);
matDVETRANS  = valDELTIME.*matVEHUVW(vecDVEVEHICLE,:);

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