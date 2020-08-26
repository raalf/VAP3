function [SURF, OUTP, COND, INPU] = fcnFORCEINTERP(SURF, OUTP, COND, INPU, VEHI, valDY, temp_y)
% Function to interpolate forces at structural nodes


[ledves, ~, ~] = find(SURF.vecDVELE > 0);
lepanels = SURF.vecDVEPANEL(ledves);

isCurWing = SURF.vecWINGTYPE(ledves) == 1;

idxdve = uint16(ledves(isCurWing));
idxpanel = lepanels(isCurWing);

m = INPU.vecM(idxpanel);

m = m(1);

% Matrix of how much we need to add to an index to get the next chordwise element
% It is done this way because n can be different for each panel. Unlike in the wake,
% we can't just add a constant value to get to the same spanwise location in the next
% row of elements
tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);

rows = repmat(idxdve,1,m) + uint16(tempm);

idx1 = SURF.vecDVELE == 1;

A(idx1) = SURF.matCOEFF(idx1,1);
B(idx1) = SURF.matCOEFF(idx1,2);
C(idx1) = SURF.matCOEFF(idx1,3);

idx2 = SURF.vecDVELE == 1; %idx2 since we need to do this even for triangles
dvenum = find(idx2==0); %dvenum in question
idxf = SURF.matADJE((ismember(SURF.matADJE(:,1), dvenum) & SURF.matADJE(:,2) == 1),3); %upstream dve num
A(idx2 ==0) = (SURF.matCOEFF(idx2==0,1)-SURF.matCOEFF(idxf,1));
B(idx2 ==0) = (SURF.matCOEFF(idx2==0,2)-SURF.matCOEFF(idxf,2));
C(idx2 ==0) = (SURF.matCOEFF(idx2==0,3)-SURF.matCOEFF(idxf,3));

% Store effective A, B, and C values for each DVE
SURF.vecDVEA = A';
SURF.vecDVEB = B';
SURF.vecDVEC = C';

etastar = ((SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1).^2).*SURF.vecDVEB(SURF.vecWINGTYPE == 1))./(3*SURF.vecDVEA(SURF.vecWINGTYPE == 1) + (SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1).^2).*SURF.vecDVEC(SURF.vecWINGTYPE == 1));
etastar_avg = mean(etastar(rows),2);
SURF.center_dist = SURF.center_dist - etastar_avg;
OUTP.vecLIFTDIST = (sum(SURF.vecDVENFREE(rows),2) + sum(SURF.vecDVENIND(rows),2)).*SURF.matLIFTDIR(isCurWing,:).*COND.valDENSITY;

chord_dist = interp1(SURF.center_dist,sum(2*SURF.vecDVEHVCRD(rows),2),SURF.span_glob(1:length(SURF.ndist{1})),'pchip');
area_dist = chord_dist.*[(SURF.span_glob(2)-SURF.span_glob(1)); (SURF.span_glob(2:length(SURF.ndist{1}))-SURF.span_glob(1:length(SURF.ndist{1})-1))];

% vecCLDIST = SURF.ndist{1}'./(0.5*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2.*INPU.vecAREA);
% j = 1;
% vecCLDIST = [];
% for i = 1:size(rows,1)
%     first = j;
%     last = j + SURF.nspnele - 1;
%     area = area_dist(last);
%     vecCLDIST = [vecCLDIST; SURF.ndist{1}(first:last)'./(0.5*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2.*area)];
%     j = j + SURF.nspnele;
% end
% 
% tempCLDIST = interp1(SURF.span_glob,vecCLDIST,temp_y,'pchip');

vecCLDIST = OUTP.vecLIFTDIST(:,3)./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2)));

tempCLDIST = interp1(SURF.center_dist,vecCLDIST,temp_y,'pchip');

chord_dist = interp1(SURF.center_dist,sum(2*SURF.vecDVEHVCRD(rows),2),temp_y,'pchip');
area_dist = chord_dist*valDY;

OUTP.vecLIFTDIST = 0.5*tempCLDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2;

tempEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),SURF.matCENTER(:,2));
tempEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matCENTER(:,2));
tempEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),SURF.matCENTER(:,2));

etastar_glob = fcnSTARGLOB([zeros(size(etastar,1),1),etastar,zeros(size(etastar,1),1)],SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1));
lemidpt = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2 - etastar_glob;

lemidpt = fcnGLOBSTAR(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1));
lemidpt(:,1) = lemidpt(:,1) + SURF.vecDVEHVCRD(SURF.vecWINGTYPE == 1)./2;
lemidpt = fcnSTARGLOB(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1));

delX = lemidpt - tempEA(1:max(max(rows)),:);

app_force = (SURF.vecDVENFREE(1:max(max(rows))) + SURF.vecDVENIND(1:max(max(rows)))).*SURF.matNORMDIR(1:max(max(rows)),:);

OUTP.vecMOMDIST_NEW = cross(delX,app_force.*COND.valDENSITY);

OUTP.vecMOMDIST_NEW = sum(reshape(OUTP.vecMOMDIST_NEW(:,2),[size(rows,1),size(rows,2)]),2);

vecCMDIST = OUTP.vecMOMDIST_NEW./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2)).*(sum(2*SURF.vecDVEHVCRD(rows),2)));
tempCMDIST = interp1(SURF.center_dist,vecCMDIST,temp_y,'pchip');

OUTP.vecMOMDIST_NEW = 0.5*tempCMDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2.*chord_dist;
end