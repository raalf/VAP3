function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITERFLEX(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

tol = 100;

iter = 1;

[OUTP] = fcnVEHENERGY(INPU, COND, SURF, OUTP, VEHI, FLAG, COND.valMAXTIME);

CL(iter,1) = OUTP.vecCL(end);
CM(iter,1) = OUTP.vecVEHCM(end);

fuseV = fcnSTARGLOB([0 0 2*OUTP.vecFUSEV(end)], 0, pi+deg2rad(COND.vecVEHPITCH),0);

CZ(iter,1) = abs((fuseV(3)+OUTP.GlobForceFuse(end,3)))/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
CX(iter,1) = abs((OUTP.GlobForce(end,1)+fuseV(1)))/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
CZtrim = (VEHI.vecFUSETAILMASS*9.81)/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);

new_tail(iter,1) = SURF.vecELEVANGLE;
new_alpha(iter,1) = COND.vecVEHALPHA;
new_fpa(iter,1) = COND.vecVEHFPA;
new_pitch(iter,1) = COND.vecVEHPITCH;

while max(abs(tol)) > 1e-3
      
    if iter > 2
        TRIM.Cmalpha = (CM(iter)-CM(iter-1))/deg2rad((new_alpha(iter)-new_alpha(iter-1)));
    end
    
    iter = iter + 1;   
            
    if iter > 2
        Mstruct = 2*OUTP.vecFUSEM(end,1) + 2*OUTP.vecFUSEV(end,1)*(OUTP.vecCGLOC(end,1) - SURF.matBEAMLOC(1,1,end));
        [OUTP] = fcnTAILMOM(SURF, VEHI, OUTP, INPU, COND, FLAG);
        M = Mstruct + OUTP.vecVEHPITCHMOM;
        vecCM = M/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC);
        new_alpha(iter,1) = COND.vecVEHALPHA + (CZtrim-CZ(iter-1,1))/((CZ(iter-1)-CZ(iter-2))/(new_alpha(iter-1)-new_alpha(iter-2)));
        new_tail(iter,1) = SURF.vecELEVANGLE + (-vecCM)/((CM(iter-1)-CM(iter-2))/(new_tail(iter-1)-new_tail(iter-2)));
    else
        new_alpha(iter,1) = COND.vecVEHALPHA + (CZtrim-CZ(iter-1,1))/(2*pi*pi/180);
        new_tail(iter,1) = SURF.vecELEVANGLE + 2;
    end
    
    % Elevator properties
    cf_c = 1; % Flap chord percentage (1 = All moving tail)
    thetaf = acos(2*cf_c - 1);
    TRIM.tau = 1 - (thetaf - sin(thetaf))/pi; % Flap effectiveness
    
    [SURF] = fcnTAILTRIM(SURF, FLAG, COND, TRIM.tau*deg2rad(new_tail(iter,1)-SURF.vecELEVANGLE), 1);

    alpha_rot = new_alpha(iter,1) - COND.vecVEHALPHA;

    COND.vecVEHALPHA = new_alpha(iter,1);
    if FLAG.GLIDING == 1
        COND.vecVEHFPA = -atand(1/(OUTP.vecCL(end)/OUTP.vecCDI(end)));
    else
        COND.vecVEHFPA = 0;
    end

    [ VEHI.matGLOBUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, COND.vecVEHALPHA, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
    pitch_rot = rad2deg(VEHI.matVEHROT(2)) - COND.vecVEHPITCH;
    COND.vecVEHPITCH = rad2deg(VEHI.matVEHROT(2)); 
    VEHI.matVEHROT(2) = deg2rad(pitch_rot);
    [SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME, SURF.matAEROCNTR] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
        INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
        SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME, SURF.matAEROCNTR);

    [ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
        SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, ~] ...
        = fcnVLST2DVEPARAM(SURF.matNPDVE, SURF.matNPVLST);

    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) - repmat(INPU.matVEHORIG(1,:),1,1);
    dcm = angle2dcm(VEHI.matVEHROT(1,3), VEHI.matVEHROT(1,1), VEHI.matVEHROT(1,2), 'ZXY');
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:)*dcm;
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) + repmat(INPU.matVEHORIG(1,:),1,1);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 0);
    
    [SURF, OUTP, INPU] = fcnBEAMFORCE(SURF, OUTP, COND, INPU, FLAG, COND.valMAXTIME);
    [OUTP] = fcnVEHENERGY(INPU, COND, SURF, OUTP, VEHI, FLAG, COND.valMAXTIME);
    
    Mstruct = 2*OUTP.vecFUSEM(end,1) + 2*OUTP.vecFUSEV(end,1)*(OUTP.vecCGLOC(end,1) - SURF.matBEAMLOC(1,1,end));
    [OUTP] = fcnTAILMOM(SURF, VEHI, OUTP, INPU, COND, FLAG);
    M = Mstruct + OUTP.vecVEHPITCHMOM;
    

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    CL(iter,1) = OUTP.vecCL(end);
    
    CM(iter,1) = M/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC);

    fuseV = fcnSTARGLOB([0 0 2*OUTP.vecFUSEV(end)], 0, pi+deg2rad(COND.vecVEHPITCH),0);

    CZ(iter,1) = abs(fuseV(3)+OUTP.GlobForceFuse(end,3))/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
    CX(iter,1) = abs(OUTP.GlobForce(end,1)+fuseV(1))/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
    
    tol = [(CZtrim - CZ(iter,1))/(CZtrim); CM(iter,1)];
    
    SURF.vecELEVANGLE = new_tail(iter);
    
    % Error handling if trim routine gets stuck for some reason
    if iter > COND.valMAXTRIMITER
        OUTP.TRIMFAIL = 1;
        return
    else
        OUTP.TRIMFAIL = 0;
    end
    
end


end