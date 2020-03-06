function [] = fcnTRIMSWEEP()

alpha_seq = [-2:1:6]';

for k = 1:length(alpha_seq)
    
    load('Discus2c_Trimmed_420CG.mat')
    
    COND.vecVEHALPHA = alpha_seq(k,1);
    
    a = alpha_seq(k,1)-alpha(end,1);

    [ VEHI.matVEHUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, a, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
    [SURF.matVLST, SURF.matCENTER, VISC.matFVLST, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNTVLST] = fcnROTVEHICLE( SURF.matDVE, SURF.matVLST, SURF.matCENTER,...
        INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, VISC.matFVLST, VISC.vecFUSEVEHICLE, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
        SURF.matNTVLST, VISC.vecFDVEVEHICLE);
    [ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
        SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, SURF.matCENTER ] ...
        = fcnVLST2DVEPARAM(SURF.matDVE, SURF.matVLST);

    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) - repmat(INPU.matVEHORIG(1,:),1,1);
    dcm = angle2dcm(VEHI.matVEHROT(1,3), VEHI.matVEHROT(1,1), VEHI.matVEHROT(1,2), 'ZXY');
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:)*dcm;
    SURF.matTRIMORIG(2,:) = SURF.matTRIMORIG(2,:) + repmat(INPU.matVEHORIG(1,:),1,1);
    
    [OUTP, SURF, COND, INPU, MISC, VEHI, VISC, FLAG, matD, vecR] = fcnVAP_MAIN(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n);
    
    Cma(k,1) = OUTP.vecVEHCM;
    CLa(k,1) = OUTP.vecCLv_AVG;
    
end

% Line of best fit to get slope
Cm_alpha = polyfit(alpha_seq,Cma,1);
Cm_alpha = Cm_alpha(1);