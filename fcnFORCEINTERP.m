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

OUTP.vecLIFTDIST = (sum(SURF.vecDVENFREE(rows),2) + sum(SURF.vecDVENIND(rows),2)).*SURF.matLIFTDIR(isCurWing,:).*COND.valDENSITY;

vecCLDIST = OUTP.vecLIFTDIST(:,3)./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2)));

tempCLDIST = interp1(SURF.center_dist,vecCLDIST,temp_y,'pchip');

chord_dist = interp1(SURF.center_dist,sum(2*SURF.vecDVEHVCRD(rows),2),temp_y,'pchip');
area_dist = chord_dist*valDY;

OUTP.vecLIFTDIST = 0.5*tempCLDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2;

tempEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),SURF.matCENTER(:,2));
tempEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matCENTER(:,2));
tempEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),SURF.matCENTER(:,2));

lemidpt = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2;

delX = lemidpt - tempEA(1:max(max(rows)),:);

lift = SURF.vecDVENFREE(1:max(max(rows))) + SURF.vecDVENIND(1:max(max(rows)));

% OUTP.vecMOMDIST_NEW = cross(delX,[zeros(max(max(rows)),2),lift(1:max(max(rows)))].*COND.valDENSITY);
OUTP.vecMOMDIST_NEW = cross(delX,lift(1:max(max(rows))).*SURF.matDVENORM(rows,:).*COND.valDENSITY);


OUTP.vecMOMDIST_NEW = sum(reshape(OUTP.vecMOMDIST_NEW(:,2),[size(rows,1),size(rows,2)]),2);

vecCMDIST = OUTP.vecMOMDIST_NEW./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2)).*(sum(2*SURF.vecDVEHVCRD(rows),2)));
tempCMDIST = interp1(SURF.center_dist,vecCMDIST,temp_y,'pchip');

OUTP.vecMOMDIST_NEW = 0.5*tempCMDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matVEHUVW).^2)))^2.*chord_dist;
end