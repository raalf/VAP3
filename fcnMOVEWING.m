function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNEWWAKEPANEL, matPANELTE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matPANELTE)
% This function moves a wing (NOT rotor) by translating all of the vertices
% in the VLST and the in-centers of each triangle in CENTER.

% INPUT:
%   valALPHA - Angle of attack for this move (radians)
%   valBETA - Sideslip angle for this move (radians)
%   valDELTIME - Timestep size for this move
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
%   matELST - List of edges from fcnTRIANG
%   vecTE - List of trailing edge edges
% OUTPUT:
%   matVLST - New vertex list with moved points
%   matCENTER - New in-center list with moved points
%   matVUINF - Freestream values at each of the vertices
%   matCUINF - Freestream values at each of the in-centers
%   matNEWWAKE - Outputs a 4 x 3 x n matrix of points for the wake DVE generation

uinf = 1;

uinf = [uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)];

translation = valDELTIME.*uinf;

% Old trailing edge vertices
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,4) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNTVLST(matDVE(vecDVETE>0,3),:);

matVLST = matVLST - repmat(translation, length(matVLST(:,1)), 1);
matCENTER = matCENTER - repmat(translation, length(matCENTER(:,1)), 1);

matNTVLST = matNTVLST - repmat(translation, length(matNTVLST(:,1)), 1);

% New trailing edge vertices
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNTVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNTVLST(matDVE(vecDVETE>0,3),:);

matNEWWAKEPANEL(:,:,1) = matPANELTE(:,:,1) - translation;
matNEWWAKEPANEL(:,:,2) = matPANELTE(:,:,2) - translation;
matNEWWAKEPANEL(:,:,3) = matPANELTE(:,:,2);
matNEWWAKEPANEL(:,:,4) = matPANELTE(:,:,1);

matPANELTE = matPANELTE - translation;

end