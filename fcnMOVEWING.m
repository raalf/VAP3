function [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecTE)
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
%   matNEWWAKE - Outputs a n x 3 x 3 matrix of points for the wake triangulation

uinf = 1;

uinf = [uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)];

translation = valDELTIME.*uinf;

% % Old trailing edge vertices
% old_te = matVLST(matDVE(vecTE,:),:);

matVLST = matVLST - repmat(translation, length(matVLST(:,1)), 1);
matCENTER = matCENTER - repmat(translation, length(matCENTER(:,1)), 1);

% % New trailing edge vertices
% new_te = matVLST(matDVE(vecTE,:),:);
% 
% % These vertices will be used to calculate the wake HDVE geometry
% matNEWWAKE(:,:,1) = [new_te(1:end/2,:); old_te((end/2)+1:end,:)];
% matNEWWAKE(:,:,2) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
% matNEWWAKE(:,:,3) = [old_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE = [];
