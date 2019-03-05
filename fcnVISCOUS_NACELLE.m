function [OUTP] = fcnVISCOUS_NACELLE(valTIMESTEP, COND, SURF, WAKE, INPU, FLAG, VISC, OUTP, VEHI)
%UNTITLED2 Summary of this function goes here
%   ONLY WORKS FOR AC FLYING IN X-DIR
% Nacelle is first computed without sym then added at the end
%%
tempLEN_DIA = 5;
% l = 1.79;
% d = 0.36;
tempNACDIAMbyPROPDIAM = 0.36/1.5240;

OUTP.valNACDIA = INPU.vecROTDIAM*tempNACDIAMbyPROPDIAM;
OUTP.valNACLEN = tempLEN_DIA*OUTP.valNACDIA;

% FF_H = 1 + 1.5/(l/d)^1.5 + 7/(l/d)^3;
% FF_T = 1 + 2.2/(l/d)^1.5 + 3.8/(l/d)^3;
% FF_JNR = 1 + 0.0025*(l/d) + 60/(l/d)^3;

% fcnINDVEL(SURF.matCENTER(idxvehrotor==1,:), valTIMESTEP, SURF, WAKE, INPU, FLAG);

tempYNACUP = INPU.matROTORHUBGLOB(:,2)+OUTP.valNACDIA/2;
tempYNACLOW = INPU.matROTORHUBGLOB(:,2)-OUTP.valNACDIA/2;
idxNACDVES = zeros(size(SURF.matCENTER,1),1);
tempNACDVES = [];
[ledves, ~, ~] = find(SURF.vecDVELE > 0);
matROWS = fcnDVEROW(ledves, SURF, INPU);

for i = 1:size(INPU.vecROTDIAM,1)
    temp = find(SURF.matCENTER(SURF.vecDVEWING==1,2)>(tempYNACLOW(i)-SURF.vecDVEHVSPN(SURF.vecDVEWING==1)) & ...
        SURF.matCENTER(SURF.vecDVEWING==1,2)<(tempYNACUP(i)+SURF.vecDVEHVSPN(SURF.vecDVEWING==1)));
	tempNACDVES = [tempNACDVES; temp];
    idxNACDVES(temp) = i;
    [matWUINF] = fcnINDVEL(SURF.matCENTER(temp,:), valTIMESTEP, SURF, WAKE, INPU, FLAG);
        
    tempCRDS = sum(2*SURF.vecDVEHVCRD(matROWS),2);
    
    tempCRDS1 = nan(size(matROWS,1),1);
    tempCRDS1(matROWS(:,1),:) = tempCRDS;
    tempCRDS = tempCRDS1;
    
    idx = temp <= size(tempCRDS,1); % NOTE: this assumes that the wing DVE is ordered first
    tempCRDS = tempCRDS(temp(idx));
    
    % Calculate taper ratio of each DVE
    idx1 = 1:(length(tempCRDS)-1); % Index of root chord elements
    idx2 = 2:length(tempCRDS); % Index of tip chord elements
    
    tempCRDS = repmat(tempCRDS,1,2);
    taper_ratio = (tempCRDS(end)./tempCRDS(1))';

    vecMAC(i,1) = (2/3)*tempCRDS(1,1).*(1 + taper_ratio + taper_ratio.^2)./(1 + taper_ratio);
    
    tempNACTOT = SURF.matUINF(temp,:) + matWUINF;
    valNACUTOT(i,1) = mean(sqrt(sum(tempNACTOT.^2,2)));
    Re(i,1) = valNACUTOT(i)*OUTP.valNACLEN(i)/VISC.valKINV;
  
% All this commented code is shit and doesnt work (each line commented represents a tear from the hard work wasted)  
%     tempVLSTL = SURF.matVLST(SURF.matDVE(temp(idx),1),:);
%     tempGLOBVLST = fcnGLOBSTAR(tempVLSTL, SURF.vecDVEROLL(temp(idx)), SURF.vecDVEPITCH(temp(idx)), SURF.vecDVEYAW(temp(idx)));
%     tempVLSTR = SURF.matVLST(SURF.matDVE(temp(idx),2),:);
%     tempGLOBVLST = [tempGLOBVLST; fcnGLOBSTAR(tempVLSTR, SURF.vecDVEROLL(temp(idx)), SURF.vecDVEPITCH(temp(idx)), SURF.vecDVEYAW(temp(idx)))];
%     
%     tempHUB = fcnGLOBSTAR(repmat(INPU.matROTORHUBGLOB(i,:), size(tempGLOBVLST,1),1), repmat(SURF.vecDVEROLL(temp(idx)),2,1), repmat(SURF.vecDVEPITCH(temp(idx)),2,1), repmat(SURF.vecDVEYAW(temp(idx)),2,1));
%     tempDEL = tempHUB-tempGLOBVLST;
%     
%     [tempVLST,IA,~] = unique([tempVLSTL(:,2);tempVLSTR(:,2)]);
%     delZ(i,1) = interp1(tempVLST,tempDEL(IA,3),tempHUB(1,2)); % Should add extrap if wing is beyond wingtip
%     
%
%     tempVLST = SURF.matVLST(SURF.matDVE(temp(idx),1:2),:);
%     [VLST,IA,~] = unique(tempVLST(:,2));
%     tempPT = INPU.matROTORHUBGLOB(i,2)-VLST;
%     idxUP = tempPT == min(tempPT(tempPT>0));
%     idxLO = tempPT == max(tempPT(tempPT<=0));
%     tempVLST = tempVLST(IA,:);
%     tempPTS = [tempVLST(idxUP,:); tempVLST(idxLO,:)];
%     
%     tempR = tempPTS(2,:)-tempPTS(1,:);
%     LEPT = [(tempR(1)*(INPU.matROTORHUBGLOB(i,2)-tempPTS(1,2))/tempR(2))+tempPTS(1,1), ...
%         INPU.matROTORHUBGLOB(i,2),...
%         (tempR(3)*(INPU.matROTORHUBGLOB(i,2)-tempPTS(1,2))/tempR(2))+tempPTS(1,3)];
%     vecL = INPU.matROTORHUBGLOB(i,:) - LEPT;
%     vecV = [SURF.matUINF(temp,1).*cosd(COND.vecVEHALPHA).*cos(COND.vecVEHBETA), ...
%         SURF.matUINF(temp,1).*sind(COND.vecVEHBETA), ...
%         SURF.matUINF(temp,1).*sind(COND.vecVEHALPHA).*cosd(COND.vecVEHBETA)];
%     vecV = vecV(IA,:);
%     vecV = vecV(idxUP,:);
%     
%     theta = acos(dot(vecV,vecL)./(norm(vecV,2).*norm(vecL,2)));
%     delZ(i,1) = sin(theta).*(norm(vecL,2));
    
%     tempSTARVLST = fcnGLOBSTAR(tempVLST, VEHI.matVEHROT(1), VEHI.matVEHROT(2), VEHI.matVEHROT(3));
%     
%     tempHUB = fcnGLOBSTAR(repmat(INPU.matROTORHUBGLOB(i,:),size(tempSTARVLST,1),1), VEHI.matVEHROT(1), VEHI.matVEHROT(2), VEHI.matVEHROT(3));
%     tempDEL = tempHUB-tempSTARVLST;
    
%     delZ(i,1) = interp1(VLST,tempDEL(IA,3),tempHUB(1,2)); % Should add extrap if wing is beyond wingtip

end

% CF = 0.074/(Re^0.2) - 1742/Re; % Turbulent + Laminar
CF = 0.074./(Re.^0.2); % Fully turbulent

Sref = pi.*OUTP.valNACDIA.^2./4;

%% Compute wetted area of the nacelle
l_small = OUTP.valNACLEN.*0.193;
R_small = OUTP.valNACDIA.*0.706/2;

l_large = OUTP.valNACLEN-l_small;
R_large = OUTP.valNACDIA./2;

a_small = (l_small.^2 + R_small.^2)./(2.*R_small);
S_small = 2.*pi.*a_small.*((R_small-a_small).*asin(l_small./a_small)+l_small);

a_large = (l_large.^2 + R_large.^2)./(2.*R_large);
S_large = 2.*pi.*a_large.*((R_large-a_large).*asin(l_large./a_large)+l_large);

Swet = S_small + S_large + (R_large.^2*pi - R_small.^2*pi);

%% Compute FF drag on nacelle
FF_S = 1 + 2.8./(tempLEN_DIA).^1.5 + 3.8./(tempLEN_DIA).^3;
CDb = CF.*FF_S.*(Swet./Sref);

OUTP.NACELLE.vecDFF = 0.5*COND.valDENSITY.*valNACUTOT.*valNACUTOT.*Sref.*CDb;

% Add high/low wing delta
vecDELCD = zeros(length(INPU.matROTORHUB(:,1)),1);
delZ = INPU.matWINGORIG(1,3)-INPU.matROTORHUB(:,3);
vecDELCD(delZ < 0) =  abs(0.008./(OUTP.valNACDIA(delZ < 0)./(delZ(delZ < 0))));
vecDELCD(delZ > 0) =  abs(0.023./(OUTP.valNACDIA(delZ > 0)./(delZ(delZ > 0))));

OUTP.NACELLE.vecHILOD = 0.5*COND.valDENSITY.*valNACUTOT.*valNACUTOT.*OUTP.valNACDIA.*vecMAC.*vecDELCD;

%% Compute interference drag
CDP = zeros(size(INPU.vecROTDIAM,1),1);
for i = 1:size(INPU.vecROTDIAM,1)
     CDP(i,1) = mean(OUTP.vecCDPDIST(idxNACDVES(ledves)==i));
end

OUTP.NACELLE.vecDINTERF = 0.5*COND.valDENSITY.*valNACUTOT.*valNACUTOT.*OUTP.valNACDIA.*vecMAC.*CDP;

%% Apply Sym

temp1 = unique([SURF.vecDVEPANEL,SURF.vecDVEWING],'rows');
vecPanelWingNumber(temp1(:,1),1) = temp1(:,2);
vecSYMWING = false(max(vecPanelWingNumber),1);
for j = 1:max(vecPanelWingNumber)
    vecSYMWING(j,1) = any(INPU.vecSYM(vecPanelWingNumber == j));
end

if vecSYMWING(1) == 1
    OUTP.NACELLE.vecDINTERF = OUTP.NACELLE.vecDINTERF*2;
    OUTP.NACELLE.vecDFF = OUTP.NACELLE.vecDFF*2;
    OUTP.NACELLE.vecHILOD = OUTP.NACELLE.vecHILOD*2;
end

% Sum all nacelle drag
OUTP.NACELLE.vecNACDRAG = OUTP.NACELLE.vecDINTERF+OUTP.NACELLE.vecDFF + OUTP.NACELLE.vecHILOD ;

