function [COND, INPU, OUTP, MISC, SURF, FLAG, TRIM, VEHI] = fcnFLEXVEHI(INPU, COND, SURF, OUTP, FLAG, MISC, VEHI, TRIM, valTIMESTEP)
% Applies gust after aeroelastic convergence
% Michael A. D. Melville, Denver, CO, 80218

% Runs stagger-step dynamic loop if valTIMESTEP > number of stiff steps,
% regardless if flight dynamics are being calculated.
if (valTIMESTEP > COND.valSTIFFSTEPS && FLAG.FLIGHTDYN == 1) || (valTIMESTEP > COND.valSTIFFSTEPS && FLAG.FLIGHTDYN == 0 && FLAG.STATICAERO == 0)
    
    % Compute structural response over stagger step range with current aero
    % forces
    COND.valSTAGGERSTEPS = floor(COND.valDELTIME/COND.valSDELTIME);
    for tempTIME = 1:COND.valSTAGGERSTEPS
        [OUTP] = fcnELASTICWING_STAGGER(OUTP, INPU, SURF, COND, VEHI, FLAG, valTIMESTEP, tempTIME);
    end
        
    if COND.valSTAGGERSTEPS > 1
        OUTP.matDEF_old = OUTP.matDEF(COND.valSTAGGERSTEPS-1:COND.valSTAGGERSTEPS,:);
        OUTP.matTWIST_old = OUTP.matTWIST(COND.valSTAGGERSTEPS-1:COND.valSTAGGERSTEPS,:);
    end
    
    OUTP.matDEF = OUTP.matDEF_old;
    OUTP.matTWIST = OUTP.matTWIST_old;
    
    OUTP.dyn_iter = OUTP.dyn_iter + 1;
    
% Runs static aeroelastic convergence loop
else
    
    tempTIME = 1;
    tol = 100;
    tol1 = 100;
    tol2 = 100;
    stability = 1;
    
    while max(abs(tol)) > 1e-7
        
        [OUTP] = fcnELASTICWING_STAGGER(OUTP, INPU, SURF, COND, VEHI, FLAG, valTIMESTEP, tempTIME);
                
        % Error checking for unstable solution
        if any(isnan(OUTP.vecDEF) == 1) || any(OUTP.vecDEF > INPU.vecSPAN)
            fprintf('\nUnstable structure solution. Reducing time step size.\n\n')
            COND.valSDELTIME = COND.valSDELTIME*0.5;
            tempTIME = 0;
            OUTP.matDEF(isnan(OUTP.matDEF)) = 0;
            OUTP.matTWIST(isnan(OUTP.matDEF)) = 0;
            OUTP.matDEF = zeros(COND.valSTIFFSTEPS,INPU.valNSELE+4);
            OUTP.matTWIST = zeros(COND.valSTIFFSTEPS,INPU.valNSELE+3);
            stability = stability + 1;
        end
        
        if stability > 10 || tempTIME > 1e5
            OUTP.TRIMFAIL = 1;
            return;
        end
        
        % Store elastic deformation results to check for convergence
        if tempTIME >= 1
            temp_def(tempTIME,1) = OUTP.matDEF(tempTIME,end-2);
            temp_twist(tempTIME,1) = OUTP.matTWIST(tempTIME,end-2);
        end
        
        if tempTIME > 5
            tol1 = (temp_def(tempTIME,1)-temp_def(tempTIME-1,1))/temp_def(tempTIME-1,1);
            tol2 = (temp_twist(tempTIME,1)-temp_twist(tempTIME-1,1))/temp_twist(tempTIME-1,1);
        end
        
        tol = [tol1; tol2];
    
        tempTIME = tempTIME + 1;
    end
    
    OUTP.matDEF = [OUTP.matDEF(end-1,:); OUTP.matDEF(end,:)];
    OUTP.matTWIST = [OUTP.matTWIST(end-1,:); OUTP.matTWIST(end,:)];
    OUTP.valSTRUCTITER = tempTIME;
    OUTP.matDEF_old = OUTP.matDEF;
    OUTP.matTWIST_old = OUTP.matTWIST;
    
    OUTP.matDEF_RELAX(OUTP.aero_iter,:) = OUTP.matDEF(end,:);
    OUTP.matTWIST_RELAX(OUTP.aero_iter,:) = OUTP.matTWIST(end,:);
    
end

% Move wing based on elastic deformations and flight dynamics. Rebuild DVE
% properties
[SURF, MISC, COND, INPU, VEHI, OUTP] = fcnMOVEFLEXVEHI2(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP);
    
[ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
    SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
    = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);

% New non-planar trailing edge vertices (used to calculate matWADJE)
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);

% Add in gust velocities to SURF.matUINF if convergence tolerance is met
if (FLAG.FLIGHTDYN == 1 && valTIMESTEP >= COND.valSTIFFSTEPS + 1)  || (FLAG.FLIGHTDYN == 0 && valTIMESTEP > COND.valSTIFFSTEPS + 1) || COND.valGUSTTIME > 1
    
%     FLAG.STATICAERO = 1;
    
    if COND.valGUSTTIME == 1
        COND.valDELTIME_old = COND.valDELTIME;
    end
    
    if COND.valGUSTTIME <= 1
        
        [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_t, SURF.matCENTER, COND.valDELTIME, valTIMESTEP);
        SURF.matUINF(:,1) = -VEHI.matGLOBUVW(1);
        [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME_old,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old,COND.start_loc,valTIMESTEP,SURF.matGUSTFIELD,SURF.vk_gust);
        OUTP.matGUSTVEL(:,valTIMESTEP) = SURF.gust_vel_old;
        COND.valGUSTTIME = COND.valGUSTTIME + 1;
        
    else
       
        % Add elastic velocities as well as gust velocity
        [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_t, SURF.matCENTER, COND.valDELTIME, valTIMESTEP);
        SURF.matUINF(:,1) = -VEHI.matGLOBUVW(1);
        [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME_old,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old,COND.start_loc,valTIMESTEP,SURF.matGUSTFIELD,SURF.vk_gust);
        OUTP.matGUSTVEL(:,valTIMESTEP) = SURF.gust_vel_old;
        COND.valGUSTTIME = COND.valGUSTTIME + 1;
       
    end
    
end

end