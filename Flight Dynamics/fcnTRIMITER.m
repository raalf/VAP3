function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

tol = 100;

iter = 1;

while max(abs(tol)) > 1e-6
    
    CL(iter,1) = OUTP.vecCL(end);
    CM(iter,1) = OUTP.vecVEHCM(end);
    
    new_tail(iter,1) = SURF.vecELEVANGLE;
    new_alpha(iter,1) = COND.vecVEHALPHA;
    
    if iter > 2
        TRIM.Cmalpha = (CM(iter)-CM(iter-1))/deg2rad((new_alpha(iter)-new_alpha(iter-1)));
    end
    
    iter = iter + 1;
        
    if iter > 2
        new_alpha(iter,1) = COND.vecVEHALPHA + (COND.CLtrim-OUTP.vecCL(end))/((CL(iter-1)-CL(iter-2))/(new_alpha(iter-1)-new_alpha(iter-2)));
        new_tail(iter,1) = SURF.vecELEVANGLE + (-OUTP.vecVEHCM(end))/((CM(iter-1)-CM(iter-2))/(new_tail(iter-1)-new_tail(iter-2)));
    else
        new_alpha(iter,1) = COND.vecVEHALPHA + (COND.CLtrim-OUTP.vecCL(end))/(2*pi*pi/180);
        new_tail(iter,1) = SURF.vecELEVANGLE + 2;
    end
    
    % Elevator properties
    cf_c = 0.25; % Flap chord percentage (1 = All moving tail)
    thetaf = acos(2*cf_c - 1);
    TRIM.tau = 1 - (thetaf - sin(thetaf))/pi; % Flap effectiveness
    
    [SURF] = fcnTAILTRIM(SURF, FLAG, COND, TRIM.tau*deg2rad(new_tail(iter,1)-SURF.vecELEVANGLE), 1);

    alpha_rot = new_alpha(iter,1) - COND.vecVEHALPHA;

    COND.vecVEHALPHA = new_alpha(iter,1);

    [ VEHI.matVEHUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, alpha_rot, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
    [SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
        INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
        SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME);

    [ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
        SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, ~] ...
        = fcnVLST2DVEPARAM(SURF.matNPDVE, SURF.matNPVLST);

    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) - repmat(INPU.matVEHORIG(1,:),1,1);
    dcm = angle2dcm(VEHI.matVEHROT(1,3), VEHI.matVEHROT(1,1), VEHI.matVEHROT(1,2), 'ZXY');
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:)*dcm;
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) + repmat(INPU.matVEHORIG(1,:),1,1);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    tol = [(COND.CLtrim - OUTP.vecCL(end))*(COND.CLtrim - OUTP.vecCL(end)); OUTP.vecVEHCM(end)*OUTP.vecVEHCM(end)];
    
    SURF.vecELEVANGLE = new_tail(iter);
    
end



end