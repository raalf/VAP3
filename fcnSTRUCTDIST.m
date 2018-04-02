function [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF)
%% Geometric Properties

% Find DVEs on LE and TE of wing
[ledves, ~, ~] = find(SURF.vecDVELE > 0);
[tedves, ~, ~] = find(SURF.vecDVETE > 0);

[matROWS] = fcnDVEROW(ledves, SURF, INPU);

tempDVEEDGECRD = abs(SURF.matNTVLST(SURF.matNTDVE(ledves,1),:) - SURF.matNTVLST(SURF.matNTDVE(tedves,4),:));
tempDVEEDGECRD = [tempDVEEDGECRD; abs(SURF.matNTVLST(SURF.matNTDVE(ledves(end),2),:) - SURF.matNTVLST(SURF.matNTDVE(tedves(end),3),:))];

tempDVEEDGECRD = sqrt(sum(tempDVEEDGECRD.^2,2)); % Chord length at each DVE edge

% Use transformation matrix to determine X,Y,Z coordinates of aerodynamic
% center based on DVE edge chord
matDVEEDGECRD = [tempDVEEDGECRD, zeros(length(tempDVEEDGECRD),2)]; % Chord vector at each DVE edge in local DVE frame (vector pointing from LE to TE)

matQTRCRD = fcnSTARGLOB(0.25*matDVEEDGECRD,[SURF.vecDVEROLL(matROWS(:,1));SURF.vecDVEROLL(matROWS(end,1))],[SURF.vecDVEPITCH(matROWS(:,1));SURF.vecDVEPITCH(matROWS(end,1))],[SURF.vecDVEYAW(matROWS(:,1));SURF.vecDVEYAW(matROWS(end,1))]);

tempLE = [SURF.matVLST(SURF.matDVE(ledves,1),:); SURF.matVLST(SURF.matDVE(ledves(end),2),:)]; % DVE LE coordinates

SURF.matAEROCNTR = tempLE + matQTRCRD; % X, Y, Z location of aerodynamic center

% Determine spanwise location (y coordinate) of DVEs
tempSPANDIST = 2*SURF.vecDVEHVSPN(ledves);
matSPANDIST = repmat(tempSPANDIST,1,size(tempSPANDIST));

tempSPAN = triu(matSPANDIST);  % Upper triangular matrix of DVE spans

SURF.vecSPANDIST = sum(tempSPAN(:,1:size(matSPANDIST,2)))';
SURF.vecSPANDIST = [0; SURF.vecSPANDIST]; % Adding location of wing root

% Calculating mean aerodynamic chord at each spanwise station to use for
% pitching moment calculations later on

idx1 = 1:(length(tempDVEEDGECRD)-1); % Index of root chord elements
idx2 = 2:length(tempDVEEDGECRD); % Index of tip chord elements

% Calculate taper ratio of each DVE
tempDVEEDGECRD = repmat(tempDVEEDGECRD,1,2);
taper_ratio = (tempDVEEDGECRD(idx2)./tempDVEEDGECRD(idx1))';

SURF.vecMAC = (2/3)*tempDVEEDGECRD(idx1,1).*(1 + taper_ratio + taper_ratio.^2)./(1 + taper_ratio); % Vector of mean aerodynamic chord at each spanwise station


%% Structural Properties

% Spanwise bending stiffness distribution. Cols 2 and 3 are the first and
% second derivatives
INPU.matEIx(:,1) = (INPU.vecEIxCOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecEIxCOEFF(2).*SURF.vecSPANDIST + INPU.vecEIxCOEFF(3))';
INPU.matEIx(:,2) = (2*INPU.vecEIxCOEFF(1).*SURF.vecSPANDIST + INPU.vecEIxCOEFF(2))';
INPU.matEIx(:,3) = (repmat(2*INPU.vecEIxCOEFF(1),1,size(SURF.vecSPANDIST)))';

% Spanwise torsional stiffness distribution. Col 2 is the first derivative
INPU.matGJt(:,1) = (INPU.vecGJtCOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecGJtCOEFF(2).*SURF.vecSPANDIST + INPU.vecGJtCOEFF(3))';
INPU.matGJt(:,2) = (2*INPU.vecGJtCOEFF(1).*SURF.vecSPANDIST + INPU.vecGJtCOEFF(2))';

INPU.vecEA = INPU.vecEACOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecEACOEFF(2).*SURF.vecSPANDIST + INPU.vecEACOEFF(3);
INPU.vecCG = INPU.vecCGCOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecCGCOEFF(2).*SURF.vecSPANDIST + INPU.vecCGCOEFF(3);
INPU.vecJT = INPU.vecJTCOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecJTCOEFF(2).*SURF.vecSPANDIST + INPU.vecJTCOEFF(3);
INPU.vecLM = INPU.vecLMCOEFF(1).*SURF.vecSPANDIST.^2 + INPU.vecLMCOEFF(2).*SURF.vecSPANDIST + INPU.vecLMCOEFF(3);

% INPU.vecLM(end) = 0; % Add zero mass to wing tip

% Determining X,Y,Z location of elastic axis (shear center) and center of
% mass (CG)
tempEA = [INPU.vecEA, zeros(length(INPU.vecEA),2)]; % Distance to EA from LE in local coordinates

tempCG = [INPU.vecCG, zeros(length(INPU.vecCG),2)]; % Distance to CG from LE in local coordinates

matCG = fcnSTARGLOB(tempCG,[SURF.vecDVEROLL(matROWS(:,1));SURF.vecDVEROLL(matROWS(end,1))],[SURF.vecDVEPITCH(matROWS(:,1));SURF.vecDVEPITCH(matROWS(end,1))],[SURF.vecDVEYAW(matROWS(:,1));SURF.vecDVEYAW(matROWS(end,1))]); % Transform to global coordinates

SURF.matSC = fcnSTARGLOB(tempEA,[SURF.vecDVEROLL(matROWS(:,1));SURF.vecDVEROLL(matROWS(end,1))],[SURF.vecDVEPITCH(matROWS(:,1));SURF.vecDVEPITCH(matROWS(end,1))],[SURF.vecDVEYAW(matROWS(:,1));SURF.vecDVEYAW(matROWS(end,1))]); % Transform to global coordinates

SURF.matSC = tempLE + SURF.matSC; % Add LE coordinates to have absolute location

matCG = tempLE + matCG;

matLSM = SURF.matSC - matCG;

matLSAC = SURF.matAEROCNTR - SURF.matSC;

SURF.vecLSAC = sqrt(sum(matLSAC.^2,2));

SURF.vecLSM = -1.*sign(matLSM(:,1)).*sqrt(sum(matLSM.^2,2)); % If +ve --> CG is ahead of SC; If -ve --> CG is behind SC

% Determine distance of each vertex to SC (to be used for twist velocity
% later)
temp_leftV = [SURF.matDVE(matROWS,1),SURF.matDVE(matROWS,4)];
temp_leftV = reshape(temp_leftV,sum(INPU.vecN,1),[]);

[move_row,~] = find(temp_leftV); % Vector corresponding to which shear center coordinate to use

tempSCLST = zeros(size(SURF.matVLST,1),3);

tempSCLST(temp_leftV,:) = SURF.matSC(move_row,:);

temp_rightV = [SURF.matDVE(matROWS,2), SURF.matDVE(matROWS,3)];
temp_rightV = reshape(temp_rightV,sum(INPU.vecN,1),[]);

[move_row,~] = find(temp_rightV); % Vector corresponding to which shear center coordinate to use
tempSCLST(temp_rightV,:) = SURF.matSC(move_row+1,:);

SURF.matSCLST = tempSCLST;

SURF.matSCLST = SURF.matSCLST - SURF.matVLST; % Matrix of vectors between shear center and vertex

figure(4)
clf
patch('Faces',SURF.matDVE,'Vertices',SURF.matVLST,'FaceColor','r')
hold on
plot3(SURF.matAEROCNTR(:,1), SURF.matAEROCNTR(:,2), SURF.matAEROCNTR(:,3),'-ok')
plot3(SURF.matSC(:,1), SURF.matSC(:,2), SURF.matSC(:,3),'-ob')
plot3(matCG(:,1), matCG(:,2), matCG(:,3),'-^g')
axis equal

end