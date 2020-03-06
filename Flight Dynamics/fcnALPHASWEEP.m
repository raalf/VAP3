function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnALPHASWEEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

alpha_input = COND.vecVEHALPHA;
alpha_seq = [COND.vecVEHALPHA; COND.vecVEHALPHA + 2];

for k = 1:length(alpha_seq)
    
    % How much to rotate the current geometry to achieve the desired AoA
    if k == 1
        alpha_rot = alpha_seq(k)-alpha_input;
    else
        alpha_rot = alpha_seq(k) - alpha_seq(k-1);
    end
    
    COND.vecVEHALPHA = alpha_seq(k);

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
    
    Cma(k,1) = OUTP.vecVEHCM;
    CLa(k,1) = OUTP.vecCL(end);
    CDa(k,1) = OUTP.vecCD(end);
    
end

% Line of best fit to get slope
CL_alpha = polyfit(alpha_seq*pi/180,CLa,1);
TRIM.CLalpha = CL_alpha(1);

TRIM.alpha0 = (-CL_alpha(2)/TRIM.CLalpha)*180/pi; % Zero-lift AoA for vehicle

CD_alpha = polyfit((alpha_seq-TRIM.alpha0)*pi/180,CDa,1);
TRIM.CDalpha = CD_alpha(1);
TRIM.valCD0 = CD_alpha(2); % Zero-lift drag coefficient

Cm_alpha = polyfit((alpha_seq-TRIM.alpha0)*pi/180,Cma,1);
TRIM.Cm0 = Cm_alpha(2);
TRIM.Cmalpha = Cm_alpha(1);

end