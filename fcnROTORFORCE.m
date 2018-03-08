function [vecCT, vecCP] = fcnROTORFORCE(matROTORDP, en, vecDVENFREE, vecDVENIND, inddrag, matUINF, vecDVEROTOR, matVEHROT, matROTORAXIS, vecROTORRPM, vecROTDIAM, matUINFROT, vecQARM, vecDVEVEHICLE, vecDELNDIST)
% Computes the thrust and power coefficients of each rotor

% Thrust direction in global reference frame
et = zeros(size(vecDVEROTOR,1),3);
tempROTORAXISDVE = zeros(size(vecDVEROTOR,1),3);
tempROTORAXISDVE(vecDVEROTOR>0,:) = matROTORAXIS(nonzeros(vecDVEROTOR),:);
for i = 1:max(vecDVEVEHICLE)
        et(vecDVEVEHICLE==i,:) = fcnSTARGLOB(tempROTORAXISDVE(vecDVEVEHICLE==i,:), matVEHROT(i,1), matVEHROT(i,2), matVEHROT(i,3));
end
% Torque direction
eq = matUINFROT./(sqrt(matUINFROT(:,1).^2+matUINFROT(:,2).^2+matUINFROT(:,3).^2));

% % Make induced drag same length as other distributions (for m>1)
% tempDi = zeros(size(vecDVENFREE,1),1);
% tempDi(vecDVETE==3) = inddrag;

% Velocitiy direction
matVELDIR = matUINF./(sqrt(matUINF(:,1).^2+matUINF(:,2).^2+matUINF(:,3).^2));

% Force distributions
vecTHRUSTDIST = dot(vecDVENFREE.*en,et,2) + dot(vecDVENIND.*en,et,2) + dot(inddrag.*matVELDIR,et,2) + dot(matROTORDP,et,2) + dot(vecDELNDIST.*en,et,2);
vecTORQUEDIST = vecQARM.*(dot(vecDVENFREE.*en,eq,2)) + vecQARM.*(dot(vecDVENIND.*en,eq,2)) + vecQARM.*(dot(inddrag.*matVELDIR,eq,2)) + vecQARM.*(dot(matROTORDP,eq,2));

for i = 1:max(vecDVEROTOR)
    idx = vecDVEROTOR == i;
    thrust(i) = sum(vecTHRUSTDIST(idx));
    torque(i) = sum(vecTORQUEDIST(idx));
end
power = torque.*2.*pi.*(vecROTORRPM'./60);

% Compute coefficients in propeller convention (not rotor convention)
vecCT = thrust'./(((vecROTORRPM/60).^2).*((vecROTDIAM).^4));
vecCQ = torque'./(((vecROTORRPM/60).^2).*((vecROTDIAM).^5));
vecCP = power'./(((vecROTORRPM/60).^3).*((vecROTDIAM).^5));

end