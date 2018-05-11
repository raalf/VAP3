function [OUTP] = fcnROTORFORCE(valTIMESTEP, matROTORDP, en, inddrag, vecDELNDIST, SURF, VEHI, INPU, COND, OUTP)
% Computes the thrust and power coefficients of each rotor

% Thrust direction in global reference frame
et = zeros(size(SURF.vecDVEROTOR,1),3);
tempROTORAXISDVE = zeros(size(SURF.vecDVEROTOR,1),3);
tempROTORAXISDVE(SURF.vecDVEROTOR>0,:) = INPU.matROTORAXIS(nonzeros(SURF.vecDVEROTOR),:);
for i = 1:max(SURF.vecDVEVEHICLE)
        et(SURF.vecDVEVEHICLE==i,:) = fcnSTARGLOB(tempROTORAXISDVE(SURF.vecDVEVEHICLE==i,:), VEHI.matVEHROT(i,1), VEHI.matVEHROT(i,2), VEHI.matVEHROT(i,3));
end
% Torque direction
eq = SURF.matUINFROT./(sqrt(SURF.matUINFROT(:,1).^2+SURF.matUINFROT(:,2).^2+SURF.matUINFROT(:,3).^2));

% % Make induced drag same length as other distributions (for m>1)
% tempDi = zeros(size(vecDVENFREE,1),1);
% tempDi(vecDVETE==3) = inddrag;

% Velocitiy direction
matVELDIR = SURF.matUINF./(sqrt(SURF.matUINF(:,1).^2+SURF.matUINF(:,2).^2+SURF.matUINF(:,3).^2));

% Force distributions
vecTHRUSTDIST = dot(SURF.vecDVENFREE.*en,et,2) + dot(SURF.vecDVENIND.*en,et,2) + dot(inddrag.*matVELDIR,et,2) + dot(matROTORDP,et,2) + dot(vecDELNDIST.*en,et,2);
vecTORQUEDIST = SURF.vecQARM.*(dot(SURF.vecDVENFREE.*en,eq,2)) + SURF.vecQARM.*(dot(SURF.vecDVENIND.*en,eq,2)) + SURF.vecQARM.*(dot(inddrag.*matVELDIR,eq,2)) + SURF.vecQARM.*(dot(matROTORDP,eq,2));
vecINVISCID_TORQUEDIST = SURF.vecQARM.*(dot(SURF.vecDVENFREE.*en,eq,2)) + SURF.vecQARM.*(dot(SURF.vecDVENIND.*en,eq,2)) + SURF.vecQARM.*(dot(inddrag.*matVELDIR,eq,2));

for i = 1:max(SURF.vecDVEROTOR)
    idx = SURF.vecDVEROTOR == i;
    thrust(i) = sum(vecTHRUSTDIST(idx));
    torque(i) = sum(vecTORQUEDIST(idx));
    inviscid_torque(i) = sum(vecINVISCID_TORQUEDIST(idx));
    
    num_blades = max(SURF.vecDVEROTORBLADE(idx));
    OUTP.ROTOR(i).vecTHRUSTDIST(valTIMESTEP,:,:) = reshape(vecTHRUSTDIST(idx), 1, [], num_blades);
    OUTP.ROTOR(i).vecTORQUEDIST(valTIMESTEP,:,:) = reshape(vecTORQUEDIST(idx), 1, [], num_blades);
end
power = torque.*2.*pi.*(COND.vecROTORRPM'./60);
inviscid_power = inviscid_torque.*2.*pi.*(COND.vecROTORRPM'./60);

% Compute coefficients in propeller convention (not rotor convention)
OUTP.vecCT(valTIMESTEP,:) = thrust'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
vecCQ = torque'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^5));
OUTP.vecCP(valTIMESTEP,:) = power'./(((COND.vecROTORRPM/60).^3).*((INPU.vecROTDIAM).^5));
OUTP.vecCPI(valTIMESTEP,:) = inviscid_power'./(((COND.vecROTORRPM/60).^3).*((INPU.vecROTDIAM).^5));

end