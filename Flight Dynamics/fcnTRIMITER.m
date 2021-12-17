function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

tol = 100;

iter = 1;

CL(iter,1) = OUTP.vecCL(end);
CM(iter,1) = OUTP.vecVEHCM(end);

CZ(iter,1) = OUTP.GlobForce(end,3)/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
CX(iter,1) = OUTP.GlobForce(end,1)/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
CZtrim = COND.vecVEHWEIGHT/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);

new_tail(iter,1) = SURF.vecELEVANGLE;
new_alpha(iter,1) = COND.vecVEHALPHA;
new_fpa(iter,1) = COND.vecVEHFPA;
new_pitch(iter,1) = COND.vecVEHPITCH;

while max(abs(tol)) > 1e-5
      
    if iter > 2
        TRIM.Cmalpha = (CM(iter)-CM(iter-1))/deg2rad((new_alpha(iter)-new_alpha(iter-1)));
    end
    
    iter = iter + 1;   
            
    if iter > 2
%         new_alpha(iter,1) = COND.vecVEHALPHA + (COND.CLtrim-OUTP.vecCL(end))/((CL(iter-1)-CL(iter-2))/(new_alpha(iter-1)-new_alpha(iter-2)));
        new_alpha(iter,1) = COND.vecVEHALPHA + (CZtrim-CZ(iter-1,1))/((CZ(iter-1)-CZ(iter-2))/(new_alpha(iter-1)-new_alpha(iter-2)));
%         new_fpa(iter,1) = COND.vecVEHFPA + (-CX(iter-1,1))/((CX(iter-1)-CX(iter-2))/(new_fpa(iter-1)-new_fpa(iter-2)));
        new_tail(iter,1) = SURF.vecELEVANGLE + (-OUTP.vecVEHCM(end))/((CM(iter-1)-CM(iter-2))/(new_tail(iter-1)-new_tail(iter-2)));
%         new_pitch(iter,1) = new_alpha(iter,1) + new_fpa(iter,1);
    else
%         new_alpha(iter,1) = COND.vecVEHALPHA + (COND.CLtrim-OUTP.vecCL(end))/(2*pi*pi/180);
        new_alpha(iter,1) = COND.vecVEHALPHA + (CZtrim-CZ(iter-1,1))/(2*pi*pi/180);
%         new_fpa(iter,1) = -atand(1/(OUTP.vecCL(end)/OUTP.vecCD(end)));
        new_tail(iter,1) = SURF.vecELEVANGLE + 2;
%         new_pitch(iter,1) = new_alpha(iter,1) + new_fpa(iter,1);
    end
    
    % Elevator properties
    cf_c = 1; % Flap chord percentage (1 = All moving tail)
    thetaf = acos(2*cf_c - 1);
    TRIM.tau = 1 - (thetaf - sin(thetaf))/pi; % Flap effectiveness
    
    [SURF] = fcnTAILTRIM(SURF, FLAG, COND, TRIM.tau*deg2rad(new_tail(iter,1)-SURF.vecELEVANGLE), 1);

    alpha_rot = new_alpha(iter,1) - COND.vecVEHALPHA;
%     pitch_rot = new_pitch(iter,1) - VEHI.vecVEHPITCH;

    COND.vecVEHALPHA = new_alpha(iter,1);
    if FLAG.GLIDING == 1
        COND.vecVEHFPA = -atand(1/(OUTP.vecCL(end)/OUTP.vecCD(end)));
%         COND.vecVEHFPA = COND.vecVEHALPHA - COND.vecVEHPITCH;
    else
        COND.vecVEHFPA = 0;
    end
%     COND.vecVEHFPA = new_fpa(iter,1);
%     VEHI.vecVEHPITCH = new_pitch(iter,1);

    [ VEHI.matGLOBUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, COND.vecVEHALPHA, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
    pitch_rot = rad2deg(VEHI.matVEHROT(2)) - COND.vecVEHPITCH;
    COND.vecVEHPITCH = rad2deg(VEHI.matVEHROT(2)); 
    VEHI.matVEHROT(2) = deg2rad(pitch_rot);
    [SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSEMASSLOC, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME, SURF.matAEROCNTR] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
        INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
        SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSEMASSLOC, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME, SURF.matAEROCNTR);

    [ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
        SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, ~] ...
        = fcnVLST2DVEPARAM(SURF.matNPDVE, SURF.matNPVLST);

    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) - repmat(INPU.matVEHORIG(1,:),1,1);
    dcm = angle2dcm(VEHI.matVEHROT(1,3), VEHI.matVEHROT(1,1), VEHI.matVEHROT(1,2), 'ZXY');
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:)*dcm;
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) + repmat(INPU.matVEHORIG(1,:),1,1);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 0);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    CL(iter,1) = OUTP.vecCL(end);
    CM(iter,1) = OUTP.vecVEHCM(end);

    CZ(iter,1) = OUTP.GlobForce(end,3)/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
    CX(iter,1) = OUTP.GlobForce(end,1)/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
    
%     tol = [(COND.CLtrim - OUTP.vecCL(end))*(COND.CLtrim - OUTP.vecCL(end)); OUTP.vecVEHCM(end)*OUTP.vecVEHCM(end)];
%     tol = [(CZtrim - CZ(iter,1))*(CZtrim - CZ(iter,1)); OUTP.vecVEHCM(end)*OUTP.vecVEHCM(end)];
    tol = [(CZtrim - CZ(iter,1))/(CZtrim); OUTP.vecVEHCM(end); CX(iter,1)];
    
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