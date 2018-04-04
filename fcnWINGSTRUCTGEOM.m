function [SURF] = fcnWINGSTRUCTGEOM(SURF, INPU)

% This function computes necessary vectors for force and moment
% distributions.
%
% INPUT:
%
% OUTPUT:
% SURF.vecSPNWSECRD - 1 x sum(INPU.vecN) vector containing the chord length at the
% mid-point of each spanwise station
% vecSPANWSEAREA - 1 x sum(INPU.vecN) vector containing the planform area at
% each spanwise station
% SURF.matQTRCRD - sum(INPU.vecN) x INPU.vecM matrix of the distance from the mid-point of
% each DVE LE to the quarter chord line

SURF.vecSPNWSECRD = [];
SURF.vecSPNWSEAREA = [];
[ledves, ~, ~] = find(SURF.vecDVELE > 0);

[matROWS] = fcnDVEROW(ledves, SURF.vecDVEPANEL, SURF.vecDVEWING, INPU.vecM, INPU.vecN);

SURF.vecSPNWSECRD = [SURF.vecSPNWSECRD; 2*sum(SURF.vecDVEHVCRD(matROWS),2)];  % Chord length at each spanwise station mid-point
SURF.vecSPNWSEAREA = [SURF.vecSPNWSEAREA; sum(SURF.vecDVEAREA(matROWS),2)];

tempMIDLE = SURF.matVLST(SURF.matDVE(:,1:2));
vecMIDLE = (tempMIDLE(:,1) + tempMIDLE(:,2))./2; % X location of mid-point on each DVE LE

matMIDLE(:,1:INPU.vecM(1)) = vecMIDLE(matROWS); % DVE LE mid-point location at each chord station

SURF.vecQTRCRD = matMIDLE(:,1) + SURF.vecSPNWSECRD*0.25; % Aerodynamic center line
SURF.matQTRCRD = -(matMIDLE - repmat(SURF.vecQTRCRD,1,INPU.vecM(1))); % Distance from DVE LE mid-point to AC

end