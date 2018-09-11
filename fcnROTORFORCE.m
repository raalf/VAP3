function [OUTP] = fcnROTORFORCE(valTIMESTEP, matROTORDP, en, inddrag, vecDELNDIST, SURF, VEHI, INPU, COND, OUTP)
% Computes the thrust and power coefficients of each rotor
TESTING = false; % TESTING SIDE FORCES AND MOMENT. WHEN FALSE CALCS ARE SKIPPED
% Thrust direction in global reference frame
et = zeros(size(SURF.vecDVEROTOR,1),3);
tempROTORAXISDVE = zeros(size(SURF.vecDVEROTOR,1),3);
tempROTORAXISDVE(SURF.vecDVEROTOR>0,:) = INPU.matROTORAXIS(nonzeros(SURF.vecDVEROTOR),:);
for i = 1:max(SURF.vecDVEVEHICLE)
        et(SURF.vecDVEVEHICLE==i,:) = fcnSTARGLOB(tempROTORAXISDVE(SURF.vecDVEVEHICLE==i,:), VEHI.matVEHROT(i,1), VEHI.matVEHROT(i,2), VEHI.matVEHROT(i,3));
end
% Torque direction
eq = SURF.matUINFROT./(sqrt(SURF.matUINFROT(:,1).^2+SURF.matUINFROT(:,2).^2+SURF.matUINFROT(:,3).^2));
if TESTING % This is under testing
    vecTHETA = (2*pi*COND.valDELTIME*COND.vecROTORRPM/60)*valTIMESTEP ...
        + (SURF.vecDVEROTORBLADE-1)*2*pi/(INPU.vecROTORBLADES);
    
    es = [-sin(vecTHETA) cos(vecTHETA) zeros(size(vecTHETA,1),1)];
    ea = [cos(vecTHETA) sin(vecTHETA) zeros(size(vecTHETA,1),1)];
end
    
% % Make induced drag same length as other distributions (for m>1)
% tempDi = zeros(size(vecDVENFREE,1),1);
% tempDi(vecDVETE==3) = inddrag;

% Velocitiy direction
matVELDIR = SURF.matUINF./(sqrt(SURF.matUINF(:,1).^2+SURF.matUINF(:,2).^2+SURF.matUINF(:,3).^2));

% Force distributions
vecTHRUSTDIST = dot(SURF.vecDVENFREE.*en,et,2) + dot(SURF.vecDVENIND.*en,et,2) + dot(inddrag.*matVELDIR,et,2) + dot(matROTORDP,et,2) + dot(vecDELNDIST.*en,et,2);
vecTORQUEDIST = SURF.vecQARM.*(dot(SURF.vecDVENFREE.*en,eq,2)) + SURF.vecQARM.*(dot(SURF.vecDVENIND.*en,eq,2)) + SURF.vecQARM.*(dot(inddrag.*matVELDIR,eq,2)) + SURF.vecQARM.*(dot(matROTORDP,eq,2));
vecINVISCID_TORQUEDIST = SURF.vecQARM.*(dot(SURF.vecDVENFREE.*en,eq,2)) + SURF.vecQARM.*(dot(SURF.vecDVENIND.*en,eq,2)) + SURF.vecQARM.*(dot(inddrag.*matVELDIR,eq,2));

if TESTING
    vecSIDEDIST = dot(SURF.vecDVENFREE.*en,es,2) + dot(SURF.vecDVENIND.*en,es,2) + dot(inddrag.*matVELDIR,es,2) + dot(matROTORDP,es,2) + dot(vecDELNDIST.*en,es,2);
    vecAXIALDIST = dot(SURF.vecDVENFREE.*en,ea,2) + dot(SURF.vecDVENIND.*en,ea,2) + dot(inddrag.*matVELDIR,ea,2) + dot(matROTORDP,ea,2) + dot(vecDELNDIST.*en,ea,2);
    vecFYDIST = vecSIDEDIST.*sin(vecTHETA) + vecAXIALDIST.*sin(pi - vecTHETA);
    vecFXDIST = vecSIDEDIST.*cos(vecTHETA) + vecAXIALDIST.*cos(pi - vecTHETA);
    vecMXDIST = vecTHRUSTDIST.*(SURF.vecQARM.*sin(vecTHETA));
    vecMYDIST = vecTHRUSTDIST.*(SURF.vecQARM.*cos(vecTHETA));
end

for i = 1:max(SURF.vecDVEROTOR)
    idx = SURF.vecDVEROTOR == i;
    thrust(i) = sum(vecTHRUSTDIST(idx));
    torque(i) = sum(vecTORQUEDIST(idx));
    inviscid_torque(i) = sum(vecINVISCID_TORQUEDIST(idx));
    
    if TESTING
        Fy(i) = sum(vecFYDIST(idx));
        Fx(i) = sum(vecFXDIST(idx));
        My(i) = sum(vecMYDIST(idx));
        Mx(i) = sum(vecMXDIST(idx));
    end
    
    num_blades = max(SURF.vecDVEROTORBLADE(idx));
    OUTP.ROTOR(i).vecTHRUSTDIST(valTIMESTEP,:,:) = reshape(vecTHRUSTDIST(idx), 1, [], num_blades);
    OUTP.ROTOR(i).vecTORQUEDIST(valTIMESTEP,:,:) = reshape(vecTORQUEDIST(idx), 1, [], num_blades);
end
power = torque.*2.*pi.*abs(COND.vecROTORRPM'./60);
inviscid_power = inviscid_torque.*2.*pi.*abs((COND.vecROTORRPM'./60));

% Compute coefficients in propeller convention (not rotor convention)
OUTP.vecCT(valTIMESTEP,:) = thrust'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
vecCQ = torque'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^5));
OUTP.vecCP(valTIMESTEP,:) = power'./((abs(COND.vecROTORRPM/60).^3).*((INPU.vecROTDIAM).^5));
OUTP.vecCPI(valTIMESTEP,:) = inviscid_power'./((abs(COND.vecROTORRPM/60).^3).*((INPU.vecROTDIAM).^5));

if TESTING
   OUTP.vecCFx(valTIMESTEP,:) = Fx'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
   OUTP.vecCFy(valTIMESTEP,:) = Fy'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
   OUTP.vecCMx(valTIMESTEP,:) = Mx'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
   OUTP.vecCMy(valTIMESTEP,:) = My'./(((COND.vecROTORRPM/60).^2).*((INPU.vecROTDIAM).^4));
end

end