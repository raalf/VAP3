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

%--------------------------------------------------------------------------

etastar = ((SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1).^2).*SURF.vecDVEB(SURF.vecWINGTYPE == 1))./(3*SURF.vecDVEA(SURF.vecWINGTYPE == 1) + (SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1).^2).*SURF.vecDVEC(SURF.vecWINGTYPE == 1));
etastar_avg = mean(etastar(rows),2); % Force location offset due to parabolic distribution

SURF.center_dist = SURF.center_dist - etastar_avg; % Adding force offset to control point y location


aero_force = (sum(SURF.vecDVENFREE(rows),2) + sum(SURF.vecDVENIND(rows),2)).*SURF.matNORMDIR(isCurWing,:).*COND.valDENSITY; % Aerodynamic force normal to beam at each DVE

% Aerodynamic force is located at the leading edge midpoint of each lifting
% line
etastar_glob = fcnSTARGLOB([zeros(size(etastar,1),1),etastar,zeros(size(etastar,1),1)],SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1)); % Force location offset due to parabolic distribution
lemidpt = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2 - etastar_glob;

lemidpt = fcnGLOBSTAR(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1));
lemidpt = fcnSTARGLOB(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == 1),SURF.vecDVEPITCH(SURF.vecWINGTYPE == 1),SURF.vecDVEYAW(SURF.vecWINGTYPE == 1));

%--------------------------------------------------------------------------

% Elastic axis location at each aerodynamic strip
% aeroEA(:,1) = interp1(SURF.matEA(:,2),SURF.matEA(:,1),SURF.matCENTER(SURF.vecWINGTYPE(ledves) == 1,2));
% aeroEA(:,2) = interp1(SURF.matEA(:,2),SURF.matEA(:,2),SURF.matCENTER(SURF.vecWINGTYPE(ledves) == 1,2));
% aeroEA(:,3) = interp1(SURF.matEA(:,2),SURF.matEA(:,3),SURF.matCENTER(SURF.vecWINGTYPE(ledves) == 1,2));
tempEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),temp_y);
tempEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),temp_y);
tempEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),temp_y);

% Elastic axis location at each structural node
% structEA(:,1) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),1),temp_y);
% structEA(:,2) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),temp_y);
% structEA(:,3) = interp1(SURF.matEALST(1:size(SURF.vecSPANLOC,1),2),SURF.matEALST(1:size(SURF.vecSPANLOC,1),3),temp_y);

% for j = 1:size(aeroEA,1)
%     if j > 1
%         idx = find(SURF.matEA(:,2) < aeroEA(j,2) & SURF.matEA(:,2) > aeroEA(j-1,2));
%         delX = aeroEA(j-1,:) - SURF.matEA(idx,:);
%         tempMOM = cross(delX,repmat(aero_force(j,:),size(delX,1),1));
%         delX = aeroEA(j,:) - SURF.matEA(idx,:);
%         tempMOM = cross(delX,repmat(aero_force(j,:),size(delX,1),1)) + tempMOM;
%         sweepMOM = [sweepMOM; tempMOM];
%     else
%         idx = find(SURF.matEA(:,2) < aeroEA(j,2));
%         delX = aeroEA(j,:) - SURF.matEA(idx,:);
%         tempMOM = cross(delX,repmat(aero_force(j,:),size(delX,1),1));
%         sweepMOM = tempMOM;
%     end
% end

chord_dist = interp1(SURF.center_dist,sum(2*SURF.vecDVEHVCRD(rows),2),SURF.span_glob(1:length(SURF.ndist{1})),'pchip');
area_dist = chord_dist.*[(SURF.span_glob(2)-SURF.span_glob(1)); (SURF.span_glob(2:length(SURF.ndist{1}))-SURF.span_glob(1:length(SURF.ndist{1})-1))];

vecCLDIST = sqrt(sum(aero_force.^2,2))./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matGLOBUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2))); % CL at each aerodynamic strip

% aero_mom_couple = cross(delX,aero_force);

tempCLDIST = interp1(SURF.center_dist,vecCLDIST,temp_y,'pchip'); % CL at each structural node

chord_dist = interp1(SURF.center_dist,sum(2*SURF.vecDVEHVCRD(rows),2),temp_y,'pchip');
area_dist = chord_dist*valDY;

OUTP.vecLIFTDIST = 0.5*tempCLDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matGLOBUVW).^2)))^2;

% delX = lemidpt - tempEA(1:max(max(rows)),:);
% 
% app_force = (SURF.vecDVENFREE(1:max(max(rows))) + SURF.vecDVENIND(1:max(max(rows)))).*SURF.matNORMDIR(1:max(max(rows)),:);
% 
% OUTP.vecMOMDIST_NEW = cross(delX,app_force.*COND.valDENSITY);
% 
% OUTP.vecMOMDIST_NEW = sum(reshape(OUTP.vecMOMDIST_NEW(:,2),[size(rows,1),size(rows,2)]),2);
% 
% vecCMDIST = OUTP.vecMOMDIST_NEW./(0.5*COND.valDENSITY*(sqrt(sum(abs(VEHI.matGLOBUVW).^2)))^2*(sum(SURF.vecDVEAREA(rows),2)).*(sum(2*SURF.vecDVEHVCRD(rows),2)));
% tempCMDIST = interp1(SURF.center_dist,vecCMDIST,temp_y,'pchip');
% 
% OUTP.vecMOMDIST_NEW = 0.5*tempCMDIST.*area_dist(1)*COND.valDENSITY*(sqrt(sum(abs(VEHI.matGLOBUVW).^2)))^2.*chord_dist;
end