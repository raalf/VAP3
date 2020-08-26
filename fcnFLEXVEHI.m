function [COND, INPU, OUTP, MISC, SURF, FLAG, TRIM, VEHI] = fcnFLEXVEHI(INPU, COND, SURF, OUTP, FLAG, MISC, VEHI, TRIM, valTIMESTEP)

% [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF);

SURF.matCENTER_old = SURF.matCENTER;

% Applies gust after aeroelastic convergence
% Michael A. D. Melville, Denver, CO, 80218
% if COND.valGUSTTIME > 1
if valTIMESTEP > COND.valSTIFFSTEPS + 1 || (valTIMESTEP >= COND.valSTIFFSTEPS && FLAG.FLIGHTDYN == 1)
    
    for tempTIME = 1:COND.valSTAGGERSTEPS
        [OUTP] = fcnELASTICWING_STAGGER2(OUTP, INPU, SURF, COND, VEHI, valTIMESTEP, tempTIME);
    end

    OUTP.matDEF_old = OUTP.matDEF(end-1:end,:);
    OUTP.matTWIST_old = OUTP.matTWIST(end-1:end,:);
    
    OUTP.dyn_iter = OUTP.dyn_iter + 1;
    
% Runs structure code until static aeroleastic convergence
else
    
    if OUTP.aero_iter == 1
        OUTP.valSTRUCTITER = COND.valSTIFFSTEPS;
    end
    
    n = 30000;
    tempTIME = 1;
    tol = 100;
    tol1 = 100;
    tol2 = 100;
    
%     OUTP.matDEF = zeros(n,INPU.valNSELE+4);
%     OUTP.matTWIST = zeros(n,INPU.valNSELE+4);
    
    while max(abs(tol)) > 1e-7
%     for tempTIME = 1:n
        
        [OUTP] = fcnELASTICWING_STAGGER(OUTP, INPU, SURF, COND, VEHI, valTIMESTEP, tempTIME);
        
        % Error checking for unstable solution
        if any(isnan(OUTP.vecDEF) == 1)
            fprintf('\nUnstable structure solution. Reducing time step size.\n\n')
            COND.valSDELTIME = COND.valSDELTIME*0.5;
            tempTIME = 0;
%             OUTP.matDEF = OUTP.matDEF*0;
%             OUTP.matTWIST = OUTP.matTWIST*0;
            OUTP.matDEF(isnan(OUTP.matDEF)) = 0;
            OUTP.matTWIST(isnan(OUTP.matDEF)) = 0;
            OUTP.matDEF = zeros(COND.valSTIFFSTEPS,INPU.valNSELE+4);
            OUTP.matTWIST = zeros(COND.valSTIFFSTEPS,INPU.valNSELE+4);
        end
        
        if tempTIME >= 1
            temp_def(tempTIME,1) = OUTP.matDEF(tempTIME,end-2);
            temp_twist(tempTIME,1) = OUTP.matTWIST(tempTIME,end-2);
        end
        
        if tempTIME > 5
            tol1 = (temp_def(tempTIME,1)-temp_def(tempTIME-1,1))/temp_def(tempTIME-1,1);
            tol2 = (temp_twist(tempTIME,1)-temp_twist(tempTIME-1,1))/temp_twist(tempTIME-1,1);
        end
        
        tol = [tol1; tol2];
    
%     end
        tempTIME = tempTIME + 1;
    end
    
    OUTP.valSTRUCTITER = tempTIME;
    OUTP.matDEF_old = OUTP.matDEF;
    OUTP.matTWIST_old = OUTP.matTWIST;
    
%     OUTP.matDEFGLOBTRIM(COND.valFULLTRIMSTEP,:) = OUTP.matDEFGLOB(end,:);
%     OUTP.matTWISTGLOBTRIM(COND.valFULLTRIMSTEP,:) = OUTP.matTWISTGLOB(end,:);
    
end

%% Relaxation of deformations

OUTP.matDEF_RELAX(OUTP.aero_iter,:) = OUTP.matDEF(end,:);
OUTP.matTWIST_RELAX(OUTP.aero_iter,:) = OUTP.matTWIST(end,:);

if OUTP.aero_iter > 1 && FLAG.FLIGHTDYN == 0
    lambda = 0.5;
    OUTP.matDEF(end-1:end,:) = repmat(lambda*(OUTP.matDEF_RELAX(OUTP.aero_iter,:)) + (1-lambda)*(OUTP.matDEF_RELAX(OUTP.aero_iter-1,:)),2,1);
    OUTP.matTWIST(end-1:end,:) = repmat(lambda*(OUTP.matTWIST_RELAX(OUTP.aero_iter,:)) + (1-lambda)*OUTP.matTWIST_RELAX(OUTP.aero_iter-1,:),2,1);

    OUTP.matDEFGLOB(valTIMESTEP,:) = lambda*(OUTP.matDEFGLOB(valTIMESTEP,:)) + (1-lambda)*OUTP.matDEFGLOB(valTIMESTEP-1,:);
    OUTP.matTWISTGLOB(valTIMESTEP,:) = lambda*(OUTP.matTWISTGLOB(valTIMESTEP,:)) + (1-lambda)*OUTP.matTWISTGLOB(valTIMESTEP-1,:);
    
    for i = 2:length(SURF.vecSPANLOC)
        OUTP.matSLOPE(valTIMESTEP,i) = asin((OUTP.matDEFGLOB(valTIMESTEP,i)-OUTP.matDEFGLOB(valTIMESTEP,i-1))/(SURF.vecSPANLOC(i)-SURF.vecSPANLOC(i-1)));
    end

end

if FLAG.FLIGHTDYN == 1 && OUTP.dyn_iter == 1
OUTP.matDEF_RELAXDYN(OUTP.dyn_iter,:) = OUTP.matDEF(end,:);
OUTP.matTWIST_RELAXDYN(OUTP.dyn_iter,:) = OUTP.matTWIST(end,:);
elseif OUTP.dyn_iter > 1 && FLAG.FLIGHTDYN == 1
    lambda = 0.5;
    OUTP.matDEF(end-1:end,:) = repmat(lambda*(OUTP.matDEF_RELAXDYN(OUTP.dyn_iter,:)) + (1-lambda)*(OUTP.matDEF_RELAXDYN(OUTP.dyn_iter-1,:)),2,1);
    OUTP.matTWIST(end-1:end,:) = repmat(lambda*(OUTP.matTWIST_RELAXDYN(OUTP.dyn_iter,:)) + (1-lambda)*OUTP.matTWIST_RELAXDYN(OUTP.dyn_iter-1,:),2,1);

    OUTP.matDEFGLOB(valTIMESTEP,:) = lambda*(OUTP.matDEFGLOB(valTIMESTEP,:)) + (1-lambda)*OUTP.matDEFGLOB(valTIMESTEP-1,:);
    OUTP.matTWISTGLOB(valTIMESTEP,:) = lambda*(OUTP.matTWISTGLOB(valTIMESTEP,:)) + (1-lambda)*OUTP.matTWISTGLOB(valTIMESTEP-1,:);
    
    for i = 2:length(SURF.vecSPANLOC)
        OUTP.matSLOPE(valTIMESTEP,i) = asin((OUTP.matDEFGLOB(valTIMESTEP,i)-OUTP.matDEFGLOB(valTIMESTEP,i-1))/(SURF.vecSPANLOC(i)-SURF.vecSPANLOC(i-1)));
    end

end

% if COND.valFULLTRIMSTEP == 1 && FLAG.FULLTRIM == 0 || FLAG.FULLTRIM == 1
    [SURF, MISC, COND, INPU, VEHI] = fcnMOVEFLEXVEHI(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP);
%     [SURF, MISC, COND] = fcnMOVEFLEXWING(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, valTIMESTEP);
% elseif COND.valFULLTRIMSTEP > 1 && FLAG.FULLTRIM == 0
%     [SURF, MISC, COND] = fcnMOVEFLEXWING(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, COND.valFULLTRIMSTEP);
% end    

[ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
    SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
    = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);

% Rotate vehicle based on flight dynamics
% if FLAG.FLIGHTDYN == 1 && valTIMESTEP > COND.valSTIFFSTEPS + 1 
%     [ ~, VEHI.matVEHROT, ~, ~] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, rad2deg(TRIM.perturb(end,4) - TRIM.perturb(end-1,4)), COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
%     [SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2,:)] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
%         INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
%         SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2,:));
% 
%     [ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
%         SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
%         = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);
% end

% New non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

if valTIMESTEP > COND.valSTIFFSTEPS + 1 || (valTIMESTEP >= COND.valSTIFFSTEPS && FLAG.FLIGHTDYN == 1)
    [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_old, SURF.matCENTER, COND.valDELTIME, COND.valSTAGGERSTEPS);
end

% Determine % relative change between aeroelastic timesteps
% tol_def = (100*abs(OUTP.matDEFGLOB(valTIMESTEP,end)-OUTP.matDEFGLOB(valTIMESTEP-COND.valSTIFFSTEPS,end))/abs(OUTP.matDEFGLOB(valTIMESTEP-COND.valSTIFFSTEPS,end)));
% tol_twist = (100*abs(OUTP.matTWISTGLOB(valTIMESTEP,end)-OUTP.matTWISTGLOB(valTIMESTEP-COND.valSTIFFSTEPS,end))/abs(OUTP.matTWISTGLOB(valTIMESTEP-COND.valSTIFFSTEPS,end)));
tol_def = 50;
tol_twist = 50;

% Add in gust velocities to SURF.matUINF if convergence tolerance is met
if (FLAG.FLIGHTDYN == 1 && valTIMESTEP > COND.valSTIFFSTEPS + 1)  || COND.valGUSTTIME > 1
    
    FLAG.STATICAERO = 1;
    
    if COND.valGUSTTIME == 1
        COND.valDELTIME_old = COND.valDELTIME;
    end
    
    if COND.valGUSTTIME <= 1
        
        [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_old, SURF.matCENTER, COND.valDELTIME, COND.valSTAGGERSTEPS);
        COND.valDELTIME = COND.valSTAGGERSTEPS*COND.valSDELTIME;
        [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME_old,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old);
        COND.valGUSTTIME = COND.valGUSTTIME + 1;
        
    else
        
        COND.valDELTIME = COND.valSTAGGERSTEPS*COND.valSDELTIME;
        
        % Add elastic velocities as well as gust velocity
        [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_old, SURF.matCENTER, COND.valDELTIME, COND.valSTAGGERSTEPS);
        [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME_old,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old);
        COND.valGUSTTIME = COND.valGUSTTIME + 1;
       
    end
    
end

end