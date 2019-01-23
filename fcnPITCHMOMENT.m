function [SURF, OUTP] = fcnPITCHMOMENT(FLAG, SURF, OUTP, INPU, COND)
% #tbt to this function. It's baaaccccccckkkkkkkkkk
% Computing the total pitching moment for the vehicle

% Find perpendicular distance from DVE control point to CG for each vehicle
for i = 1:INPU.valVEHICLES
    
    % Find LE midpoint location in xyz
    temp = fcnGLOBSTAR(SURF.matCENTER, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW);
    SURF.matDVELEMP = fcnSTARGLOB([temp(:,1) - SURF.vecDVEHVCRD,temp(:,2),temp(:,3)], SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW);
    
    % Compute moment arm from CG to DVE LE midpoint
    MomArm(:,:,i) = (INPU.vecVEHCG(i,1) - SURF.matDVELEMP(:,1)).*cos(COND.vecVEHALPHA(i))...
        + (INPU.vecVEHCG(i,3) - SURF.matCENTER(:,3)).*sin(COND.vecVEHALPHA(i));
    
    % Compute moment from lifting lines
    deltaM = (SURF.vecDVENFREE + SURF.vecDVENIND).*MomArm;
    
    % Total vehicle moment (Moment/density)
    OUTP.vecVEHMOM(i) = sum(deltaM,1) + sum(OUTP.vecCMDIST.*(0.5*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA(i)*INPU.vecCMAC(i)),1);
    
    % Non-dimensionalize
    OUTP.vecVEHCM(i) = OUTP.vecVEHMOM(i)/(0.5*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA(i)*INPU.vecCMAC(i));

end