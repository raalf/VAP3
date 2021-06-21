function [INPU, SURF] = fcnVEHISTRUCT(COND, INPU, SURF, FLAG)
%% Geometric Properties

% Find indices for flexible wing(s)
SURF.idxFLEX = find(SURF.vecDVEWING == find(FLAG.vecFLEXIBLE == 1));
SURF.idxTAIL = find(SURF.vecDVEWING == 2);

% Find DVEs on LE and TE of flexible wing
[ledves, ~, ~] = find(SURF.vecDVELE(SURF.idxFLEX) > 0);
[tedves, ~, ~] = find(SURF.vecDVETE(SURF.idxFLEX) > 0);

lepanels = SURF.vecDVEPANEL(ledves);

% Determine DVEs in each spanwise station
for i = 1:length(find(FLAG.vecFLEXIBLE == 1))

	idxdve = ledves(SURF.vecDVEWING(ledves) == i);
	idxpanel = lepanels(SURF.vecDVEWING(ledves) == i);

    m = INPU.vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
    
    matROWS{i} = repmat(idxdve,1,m) + double(tempm);

end

SURF.vecCRDDIST = sum(2.*SURF.vecDVEHVCRD(matROWS{1}),2); % Chord dist. across span at each dve ctrl pt

% Determine spanwise location (y coordinate) of DVEs
SURF.vecSPANDIST = SURF.matCENTER(matROWS{1}(:,1),2);

SURF.vecSPANLOC = [0; SURF.matNPVLST(SURF.matNPDVE(matROWS{1}(:,1),2),2)];

SURF.vecSTRUCTSPNDIST = SURF.vecSPANLOC;
INPU.valDY = SURF.vecSPANLOC(2:end)-SURF.vecSPANLOC(1:end-1);
INPU.valNSELE = length(SURF.vecSTRUCTSPNDIST);

SURF.idxNPVLST_tail = unique(SURF.matNPDVE(SURF.vecWINGTYPE == 2,:)); % Rows in VLST corresponding to tail coordinates

idx_vlst = (1:size(SURF.matNPVLST,1))~=SURF.idxNPVLST_tail;
idx_vlst = any(idx_vlst == 0,1);
idx_vlst = find(idx_vlst == 0);

SURF.idx_struct = [];
j = 1;

% Find indexing vector for matNPVLST. Columns idx_struct corresponds what index in
% the structural deflection vector should be used to move the corresponding
% indices in matNPVLST
for i = 1:length(SURF.vecSPANLOC)
    temp = SURF.matNPVLST(idx_vlst,2) - SURF.vecSPANLOC(i);
    temp2 = find(temp == 0);
    
    if length(temp2) > size(SURF.idx_struct,1) && i > 1
        SURF.idx_struct(length(temp2),:) = nan;
    elseif length(temp2) < size(SURF.idx_struct,1) && i > 1
        SURF.idx_struct(length(temp2)+1:size(SURF.idx_struct,1),end+1) = nan;
    end
    SURF.idx_struct(1:length(temp2),j) = temp2;
    j = j + 1;
end

for i = 1:length(SURF.vecSPANLOC)
    idx = find(isnan(SURF.idx_struct(:,i)) | SURF.idx_struct(:,i) == 0);

    if ~isempty(idx)
        SURF.idx_struct(idx,i) = SURF.idx_struct(idx(1)-1,i);
    end
end

struct_edgecrd = interp1(SURF.vecSPANDIST,SURF.vecCRDDIST,SURF.vecSTRUCTSPNDIST,'linear','extrap'); % Chord dist. at structural nodes

% Use transformation matrix to determine X,Y,Z coordinates of aerodynamic
% center based on DVE edge chord
matCRDDIST = [struct_edgecrd, zeros(length(struct_edgecrd),2)]; % Chord vector at each DVE edge in local DVE frame (vector pointing from LE to TE)

SURF.matQTRCRD = fcnSTARGLOB(0.25*matCRDDIST,interp1(SURF.vecSPANDIST,SURF.vecDVEROLL(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEPITCH(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEYAW(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'));

% DVE LE midpt coordinates
tempLE = [SURF.matVLST(SURF.matDVE(ledves,1),:); SURF.matVLST(SURF.matDVE(ledves(end),2),:)]; 
tempLE = (tempLE(1:end-1,:) + tempLE(2:end,:))./2;

SURF.matAEROCNTR = interp1(SURF.vecSPANDIST,tempLE,SURF.vecSTRUCTSPNDIST,'linear','extrap') + SURF.matQTRCRD; % X, Y, Z location of aerodynamic center at structure nodes

% Calculating mean aerodynamic chord at each spanwise station to use for
% pitching moment calculations later on

tempDVEEDGECRD = abs(SURF.matNPVLST(SURF.matNPDVE(ledves,1),:) - SURF.matNPVLST(SURF.matNPDVE(tedves,4),:));
tempDVEEDGECRD = [tempDVEEDGECRD; abs(SURF.matNPVLST(SURF.matNPDVE(ledves(end),2),:) - SURF.matNPVLST(SURF.matNPDVE(tedves(end),3),:))];

tempDVEEDGECRD = sqrt(sum(tempDVEEDGECRD.^2,2)); % Chord length at each DVE edge

idx1 = 1:(length(tempDVEEDGECRD)-1); % Index of root chord elements
idx2 = 2:length(tempDVEEDGECRD); % Index of tip chord elements

% Calculate taper ratio of each DVE
tempDVEEDGECRD = repmat(tempDVEEDGECRD,1,2);
taper_ratio = (tempDVEEDGECRD(idx2)./tempDVEEDGECRD(idx1))';

SURF.vecMAC = (2/3)*tempDVEEDGECRD(idx1,1).*(1 + taper_ratio + taper_ratio.^2)./(1 + taper_ratio); % Vector of mean aerodynamic chord at each spanwise station

%% Structural Properties

% Spanwise bending stiffness distribution. Cols 2 and 3 are the first and
% second derivatives
INPU.matEIx(:,1) = (INPU.vecEIxCOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecEIxCOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecEIxCOEFF(3))';

% First derivative
INPU.matEIx(2:end-1,2) = (INPU.matEIx(3:end,1)-INPU.matEIx(1:end-2,1))./(2*INPU.valDY(2:end));
INPU.matEIx(1,2) = (-3*INPU.matEIx(1,1) + 4*INPU.matEIx(2,1) - INPU.matEIx(3,1))./(2*INPU.valDY(1));
INPU.matEIx(end,2) = (3*INPU.matEIx(end,1) - 4*INPU.matEIx(end-1,1) + INPU.matEIx(end-2,1))./(2*INPU.valDY(end));
% INPU.matEIx(:,2) = (2*INPU.vecEIxCOEFF(1).*SURF.vecSTRUCTSPNDIST + INPU.vecEIxCOEFF(2))';
% INPU.matEIx(:,3) = (repmat(2*INPU.vecEIxCOEFF(1),1,size(SURF.vecSTRUCTSPNDIST)))';

% Second derivative
INPU.matEIx(2:end-1,3) = (INPU.matEIx(1:end-2,1) - 2*INPU.matEIx(2:end-1,1) + INPU.matEIx(3:end,1))./(INPU.valDY(2:end).*INPU.valDY(2:end));
INPU.matEIx(1,3) = (2*INPU.matEIx(1,1) - 5*INPU.matEIx(2,1) + 4*INPU.matEIx(3,1) - INPU.matEIx(4,1))./(INPU.valDY(1)*INPU.valDY(1)*INPU.valDY(1));
INPU.matEIx(end,3) = (2*INPU.matEIx(end,1) - 5*INPU.matEIx(end-1,1) + 4*INPU.matEIx(end-2,1) - INPU.matEIx(end-3,1))./(INPU.valDY(end)*INPU.valDY(end)*INPU.valDY(end));

% Spanwise torsional stiffness distribution. Col 2 is the first derivative
INPU.matGJt(:,1) = (INPU.vecGJtCOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecGJtCOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecGJtCOEFF(3))';
% INPU.matGJt(:,2) = (2*INPU.vecGJtCOEFF(1).*SURF.vecSTRUCTSPNDIST + INPU.vecGJtCOEFF(2))';

% First derivative
INPU.matGJt(2:end-1,2) = (INPU.matGJt(3:end,1)-INPU.matGJt(1:end-2,1))./(2*INPU.valDY(2:end));
INPU.matGJt(1,2) = (-3*INPU.matGJt(1,1) + 4*INPU.matGJt(2,1) - INPU.matGJt(3,1))./(2*INPU.valDY(1));
INPU.matGJt(end,2) = (3*INPU.matGJt(end,1) - 4*INPU.matGJt(end-1,1) + INPU.matGJt(end-2,1))./(2*INPU.valDY(end));

INPU.vecEA = INPU.vecEACOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecEACOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecEACOEFF(3);
INPU.vecCG = INPU.vecCGCOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecCGCOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecCGCOEFF(3);
INPU.vecJT = INPU.vecJTCOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecJTCOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecJTCOEFF(3);
INPU.vecLM = INPU.vecLMCOEFF(1).*SURF.vecSTRUCTSPNDIST.^2 + INPU.vecLMCOEFF(2).*SURF.vecSTRUCTSPNDIST + INPU.vecLMCOEFF(3);

% INPU.vecLM(end) = 0; % Add zero mass to wing tip

% Determining X,Y,Z location of elastic axis (shear center) and center of
% mass (CG)
tempEA = [INPU.vecEA.*struct_edgecrd(:,1), zeros(length(INPU.vecEA),2)]; % Distance to EA from LE in local coordinates

tempCG = [INPU.vecCG.*struct_edgecrd(:,1), zeros(length(INPU.vecCG),2)]; % Distance to CG from LE in local coordinates

matCG = fcnSTARGLOB(tempCG,interp1(SURF.vecSPANDIST,SURF.vecDVEROLL(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEPITCH(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEYAW(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap')); % Transform to global coordinates

SURF.matEA = fcnSTARGLOB(tempEA,interp1(SURF.vecSPANDIST,SURF.vecDVEROLL(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEPITCH(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEYAW(matROWS{1}(:,1)),SURF.vecSTRUCTSPNDIST,'linear','extrap')); % Transform to global coordinates

SURF.matEA = interp1(SURF.vecSPANDIST,tempLE,SURF.vecSTRUCTSPNDIST,'linear','extrap') + SURF.matEA; % Add LE coordinates to have absolute location

SURF.matCG = interp1(SURF.vecSPANDIST,tempLE,SURF.vecSTRUCTSPNDIST,'linear','extrap') + matCG;

matLSM = SURF.matEA - SURF.matCG;

matLSAC = SURF.matAEROCNTR - SURF.matEA;

SURF.vecLSAC = sqrt(sum(matLSAC.^2,2));

SURF.vecLSM = -1.*sign(matLSM(:,1)).*sqrt(sum(matLSM.^2,2)); % If +ve --> EA is ahead of CG; If -ve --> EA is behind CG

% SURF.valTBOOM = abs(SURF.matAEROCNTR(1,1) - SURF.matTRIMORIG(2,1)); % Tail boom length from wing elastic axis to tail qtr crd

% Compute wing mass and CG location in aeordynamic coordinate system
SURF.vecVEHMASS = interp1(SURF.vecSTRUCTSPNDIST,INPU.vecLM,SURF.matCENTER(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1,2)).*(2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)); % Wing mass values at each DVE control point y-location

SURF.vecWINGCG = interp1(SURF.vecSTRUCTSPNDIST,SURF.matCG,SURF.matCENTER(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1,2)); % Wing CG location at each DVE control point y-location


% -------------------------------------------------------------------------
temp_vecEA = INPU.vecEACOEFF(1).*SURF.vecSPANLOC.^2 + INPU.vecEACOEFF(2).*SURF.vecSPANLOC + INPU.vecEACOEFF(3);
temp_vecCG = INPU.vecCGCOEFF(1).*SURF.vecSPANLOC.^2 + INPU.vecCGCOEFF(2).*SURF.vecSPANLOC + INPU.vecCGCOEFF(3);

tempEA = [temp_vecEA.*tempDVEEDGECRD(:,1), zeros(length(temp_vecEA),2)]; % Distance to EA from LE in local coordinates
tempCG = [temp_vecCG.*tempDVEEDGECRD(:,1), zeros(length(temp_vecCG),2)]; % Distance to EA from LE in local coordinates

temp_matEA = fcnSTARGLOB(tempEA,interp1(SURF.vecSPANDIST,SURF.vecDVEROLL(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEPITCH(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEYAW(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap')); % Transform to global coordinates

temp_matCG = fcnSTARGLOB(tempCG,interp1(SURF.vecSPANDIST,SURF.vecDVEROLL(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEPITCH(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap'),...
    interp1(SURF.vecSPANDIST,SURF.vecDVEYAW(matROWS{1}(:,1)),SURF.vecSPANLOC,'linear','extrap')); % Transform to global coordinates

tempLE = [SURF.matVLST(SURF.matDVE(ledves,1),:); SURF.matVLST(SURF.matDVE(ledves(end),2),:)]; 

temp_matEA = tempLE + temp_matEA; % Add LE coordinates to have absolute location
temp_matCG = tempLE + temp_matCG; % Add LE coordinates to have absolute location

% Determine distance of each vertex to SC (to be used for twist velocity
% later)
temp_leftV = [SURF.matDVE(matROWS{1},1),SURF.matDVE(matROWS{1},4)];
temp_leftV = reshape(temp_leftV,sum(INPU.vecN(1:max(SURF.vecDVEPANEL(SURF.vecWINGTYPE == 1))),1),[]);

[move_row,~] = find(temp_leftV); % Vector corresponding to which shear center coordinate to use

tempSCLST = zeros(size(SURF.matVLST,1),3);

tempSCLST(temp_leftV,:) = temp_matEA(move_row,:);

temp_rightV = [SURF.matDVE(matROWS{1},2), SURF.matDVE(matROWS{1},3)];
temp_rightV = reshape(temp_rightV,sum(INPU.vecN(1:max(SURF.vecDVEPANEL(SURF.vecWINGTYPE == 1))),1),[]);

[move_row,~] = find(temp_rightV); % Vector corresponding to which shear center coordinate to use
tempSCLST(temp_rightV,:) = temp_matEA(move_row+1,:);

SURF.matEALST = tempSCLST;

tempCGLST = zeros(size(SURF.matVLST,1),3);

tempCGLST(temp_leftV,:) = temp_matCG(move_row,:);

temp_rightV = [SURF.matDVE(matROWS{1},2), SURF.matDVE(matROWS{1},3)];
temp_rightV = reshape(temp_rightV,sum(INPU.vecN(1:max(SURF.vecDVEPANEL(SURF.vecWINGTYPE == 1))),1),[]);

[move_row,~] = find(temp_rightV); % Vector corresponding to which shear center coordinate to use
tempCGLST(temp_rightV,:) = temp_matCG(move_row+1,:);

SURF.matCGLST = tempCGLST;

end