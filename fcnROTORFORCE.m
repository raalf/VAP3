function [ vecCT ] = fcnROTORFORCE( en, vecDVENFREE, vecDVENIND, inddrag, matUINF, vecDVEROTOR, matVEHROT, matROTORAXIS, vecROTORRPM, vecROTDIAM, matCENTER, valDELTIME, valTIMESTEP, vecCTCONV)

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

vecCT = thrust'./(((vecROTORRPM/60).^2).*((vecROTDIAM).^4));

end