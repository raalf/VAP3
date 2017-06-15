function [valNELE, matNEWNPVLST, vecAIRFOIL, vecDVELE, vecDVETE, ...
    vecDVEYAW, vecDVEPANEL, vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, ...
    vecDVEROLL, vecDVEAREA, vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, vecDVELESWP, ...
    vecDVEHVCRD, vecDVEHVSPN, vecSYM, vecQARM, matADJE, matNEWCENTER, matNEWVLST, matDVE, matNEWDVENORM, matVLST] = fcnDVEMULTIROTOR3(...
    valNELE, valNUMB, vecDVETIP, vecDVETESWP, vecDVEPITCH, vecDVEWING, vecDVEMCSWP, vecM, vecN, vecDVEPANEL, vecDVEROLL, vecDVELESWP, ...
    vecDVEYAW, vecDVEHVCRD, vecDVEHVSPN, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, vecROTAX, vecAIRFOIL, NPP, matDVE, matADJE, ...
    P, matCENTER, matDVENORM)
% This function modifies the created DVEs and all required input values
%	for multiple rotor blades.
%
% Inputs:
%   Parameters from input file and fcnGENERATEDVES that require
%   modification
%
% Outputs:
%   Each input has been modified for the given number of blades:
%   From input file
%       valAREA - multiply by valNUMB
%       vecSYM - Add (will be deleted)
%       vecAIRFOIL - add
%       vecN - add
%       venM - add
%
%   From Generate DVEs
%		matCENTER - rotate values center values
%		vecDVEHVSPN - add to vector for each blade
%		vecDVEHVCRD - add to vector for each blade
%		vecDVELESWP - add to vector for each blade
%		vecDVEMCSWP - add to vector for each blade
%		vecDVETESWP - add to vector for each blade
%		vecDVEROLL - add to vector for each blade
%		vecDVEPITCH - add to vector for each blade
%		vecDVEYAW - add to vector for each blade
%		vecDVEAREA - add to vector for each blade
%		matDVENORM - rotate values normal vectors
%		matVLST - rotate vertices
%		valNELE - multiply by valNUMB
%		matDVE - matDVE + max(matDVE) for each blade
%		matADJE - matADJE + max(matADJE) for each blade
%		vecDVESYM - delete in later versions
%		vecDVETIP - add to vector for each blade
%		vecDVEWING - vecDVEWING + max(vecDVEWING) for each blade
%		vecDVELE - add to vector for each blade
%		vecDVETE - add to vector for each blade
%		vecDVEPANEL - vecDVEPANEL + max(vecDVEPANEL) for each blade


%% Create rotation matrix
% Angle of each blade from original (+ve CCW & in rad)
tempTHETA = (0:2*pi/valNUMB:2*pi)';
tempTHETA(valNUMB+1) = [];

% Develop 
tempTOP = [cos(tempTHETA), -sin(tempTHETA), zeros([valNUMB,1])];
tempMID = [sin(tempTHETA), cos(tempTHETA), zeros([valNUMB,1])];
tempLOW = [zeros([valNUMB, 1]), zeros([valNUMB, 1]), ones([valNUMB, 1])];
tempROTATE = reshape(tempTOP',[1, 3, valNUMB]);
tempROTATE(2,:,:) = reshape(tempMID',[1, 3, valNUMB]);
tempROTATE(3,:,:) = reshape(tempLOW',[1, 3, valNUMB]);

% 2D rotation matrix
tempROTATE2D = (reshape(permute(tempROTATE,[2,1,3]),[3 valNUMB*3]))';


%% Parameters the must Rotate
% Make Each point relative to rotation axis
tempCENTER = matCENTER - vecROTAX;
P = P - vecROTAX;

% Rotate values

tempP1 = [];
tempP2 = [];
tempP3 = [];
tempP4 = [];

for i = 2:valNUMB
    tempP1 = cat(1,tempP1,(tempROTATE(:,:,i)*P(:,:,1)')');
    tempP2 = cat(1,tempP2,(tempROTATE(:,:,i)*P(:,:,2)')');
    tempP3 = cat(1,tempP3,(tempROTATE(:,:,i)*P(:,:,3)')');
    tempP4 = cat(1,tempP4,(tempROTATE(:,:,i)*P(:,:,4)')');
end

matNEWCENTER = (tempP1 + tempP2 + tempP3 + tempP4)./4;

tempNPVLST = matNPVLST - vecROTAX;
tempNEWNPVLST = (tempROTATE2D*tempNPVLST')';
temp = reshape(tempNEWNPVLST,[numel(tempNPVLST)/3,3,valNUMB]);
matNEWNPVLST = reshape(permute(temp,[2,1,3]),[3,valNUMB*(numel(tempNPVLST)/3)])' + vecROTAX;

%% Parameters to increase vector for number of blades
vecNEWAIRFOIL = repmat(vecAIRFOIL,[valNUMB-1,1]);
vecNEWSYM = repmat(vecSYM,[valNUMB-1,1]); % TO BE DELETE IN LATER VERSIONS
vecNEWN = repmat(vecN,[valNUMB-1,1]);
vecNEWM = repmat(vecM,[valNUMB-1,1]);

[ vecNEWDVEHVSPN, vecNEWDVEHVCRD, ...
    vecNEWDVEROLL, vecNEWDVEPITCH, vecNEWDVEYAW,...
    vecNEWDVELESWP, vecNEWDVEMCSWP, vecNEWDVETESWP, ...
    vecNEWDVEAREA, matNEWDVENORM, ...
    matNEWVLST, matNEWDVE, ~, idxVLST] = fcnDVECORNER2PARAM( matNEWCENTER, tempP1, tempP2, tempP3, tempP4);

% Chordwise radial distances
vecQARM = abs(matNEWCENTER(:,2)-vecROTAX(2));

fcnPLOTBODY(1, size(matNEWCENTER,1), matNEWDVE, matNEWVLST, matNEWCENTER)
end

