function [matVLST, matCENTER, vecROTORIG] = fcnMOVEROTOR(vecROTORIG, valALPHAR, valAZNUM, valJ, matVLST, matCENTER)

% Calculating the farthest point from the rotation axis, which will be valD (diameter)
valD = max(sqrt((matVLST(:,1) + vecROTORIG(1)).^2 + (matVLST(:,2) + vecROTORIG(2)).^2));

% Finding the rotation angle for the timestep
valDELROT = (2*pi)/valAZNUM;

% Finding how much to translate the rotor based on alpha
vecL = [-((valJ*valD)/valAZNUM)*cos(valALPHAR) 0 ((valJ*valD)/valAZNUM)*sin(valALPHAR)];


tempVLST = matVLST - repmat(vecROTORIG, length(matVLST(:,1)),1);
tempCENTER = matCENTER - repmat(vecROTORIG, length(matCENTER(:,1)),1);

ROT = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

vlst2 = (ROT*tempVLST')' + repmat(vecROTORIG, length(matVLST(:,1)),1) + vecL;
center2 = (ROT*tempCENTER')' + repmat(vecROTORIG, length(matCENTER(:,1)),1) + vecL;

% NEED THE TRAILING EDGE VECTOR!!!
% % Old trailing edge of wing
% old_te = matVLST(matELST(vecTE,:),:);
% 
% % New trailing edge of wing
% new_te = vlst2(matELST(vecTE,:),:);
% 
% % These vertices will be used to calculate the wake HDVE geometry
% matNEWWAKE(:,:,1) = [new_te(1:end/2,:); old_te((end/2)+1:end,:)];
% matNEWWAKE(:,:,2) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
% matNEWWAKE(:,:,3) = [old_te(1:end/2,:); new_te((end/2)+1:end,:)];

% Replacing old values with new locations
matVLST = vlst2;
matCENTER = center2;

vecROTORIG = vecROTORIG + vecL;

end

