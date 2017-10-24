function [ vecCT, vecCTCONV ] = fcnROTORFORCE( en, vecDVENFREE, vecDVENIND, inddrag, matUINF, vecDVEROTOR, matVEHROT, matROTORAXIS, vecROTORRPM, vecROTDIAM, matCENTER, valDELTIME, valTIMESTEP, vecCTCONV)

% Thrust direction in global reference frame
et = repmat(fcnSTARGLOB(matROTORAXIS, matVEHROT(:,1), matVEHROT(:,2), matVEHROT(:,3)),size(vecDVENFREE,1)/max(vecDVEROTOR),1);

% % Make induced drag same length as other distributions (for m>1)
% tempDi = zeros(size(vecDVENFREE,1),1);
% tempDi(vecDVETE==3) = inddrag;

% Velocitiy direction
matVELDIR = matUINF./(sqrt(matUINF(:,1).^2+matUINF(:,2).^2+matUINF(:,3).^2));

% Force distributions
vecTHRUSTDIST = dot(vecDVENFREE.*en,et,2) + dot(vecDVENIND.*en,et,2) + dot(inddrag.*matVELDIR,et,2);

for i = 1:max(vecDVEROTOR)
    idx = vecDVEROTOR == i;
    thrust(i) = sum(vecTHRUSTDIST(idx));
end

tempCT = thrust'./(((vecROTORRPM/60).^2).*((vecROTDIAM).^4));
vecAZNUM = (1./(abs(vecROTORRPM)/60))./valDELTIME;
if vecAZNUM < 1
    disp('Timestep size too great, error in fcnROTORFORCE.')
end
temp = valTIMESTEP - ( floor( (valTIMESTEP-1)/vecAZNUM))*(vecAZNUM);

vecCTCONV(temp,:) = tempCT';

vecCT = mean(vecCTCONV,1);

% hold on
% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),en(:,1),en(:,2),en(:,3))
% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matVELDIR(:,1),matVELDIR(:,2),matVELDIR(:,3))
% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),et(:,1),et(:,2),et(:,3))
end