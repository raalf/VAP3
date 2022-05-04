function [SURF, INPU, COND, MISC, VISC, OUTP, FLAG, TRIM, VEHI, n] = fcnMOVESTRUCTURE(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, TRIM, valTIMESTEP, n)

if (valTIMESTEP <= COND.valSTIFFSTEPS && FLAG.FLIGHTDYN == 0) || (FLAG.STIFFWING == 1 && FLAG.FLIGHTDYN == 0) || (FLAG.FLIGHTDYN == 1 && valTIMESTEP <= COND.valSTIFFSTEPS)
    
    
    if valTIMESTEP > 1
        OUTP.matDEFGLOB(valTIMESTEP,:) = OUTP.matDEFGLOB(valTIMESTEP-1,:);

        OUTP.matTWISTGLOB(valTIMESTEP,:) = OUTP.matTWISTGLOB(valTIMESTEP-1,:);
    else
        if sum(any(abs(OUTP.matDEFGLOB) > 0)) > 0
            OUTP.matDEFGLOB(valTIMESTEP,:) = OUTP.matDEFGLOB(end,:);
            OUTP.matTWISTGLOB(valTIMESTEP,:) = OUTP.matTWISTGLOB(end,:);
        else
            OUTP.matDEFGLOB(valTIMESTEP,:) = zeros(1,length(find(SURF.vecDVELE(SURF.vecWINGTYPE == 1) == 1)) + 1);
            OUTP.matTWISTGLOB(valTIMESTEP,:) = zeros(1,length(find(SURF.vecDVELE(SURF.vecWINGTYPE == 1) == 1)) + 1);
        end
    end
        
    if valTIMESTEP < 5 % Build up wake before allowing flexible wing to move
        [SURF, INPU, MISC, VISC, OUTP, VEHI] = fcnSTIFFWING(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valTIMESTEP);
    else
        [SURF, MISC, COND, INPU, VEHI, OUTP] = fcnMOVEFLEXVEHI2(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP);
        
        [ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
            SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
            = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);

        % New non-planar trailing edge vertices (used to calculate matWADJE)
        MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,1) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,4),:);
        MISC.matNPNEWWAKE(1:length(find(SURF.vecDVETE(SURF.idxFLEX) == 3)),:,2) = SURF.matNPVLST(SURF.matNPDVE(SURF.vecDVETE(SURF.idxFLEX)>0,3),:);
        
        [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_t, SURF.matCENTER, COND.valDELTIME, valTIMESTEP);
        
        if FLAG.GUSTMODE > 0
            [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old,COND.start_loc);
        end
    
    end

    
elseif FLAG.FLIGHTDYN == 1 && valTIMESTEP > COND.valSTIFFSTEPS && FLAG.STIFFWING == 1
    
    SURF.matCENTER_old = SURF.matCENTER;

    [SURF, MISC, COND, INPU, VEHI] = fcnMOVESURFACE_FLIGHTDYN(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM, valTIMESTEP);
    
    [ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
    SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, SURF.matVLST, SURF.matDVE, SURF.matCENTER, MISC.matNEWWAKE ] ...
    = fcnVLST2DVEPARAM_NEW(SURF.matNPDVE, SURF.matNPVLST, MISC.matNEWWAKE, SURF.vecDVETE);

    OUTP.matDEFGLOB(valTIMESTEP,:) = OUTP.matDEFGLOB(valTIMESTEP-1,:);
    OUTP.matTWISTGLOB(valTIMESTEP,:) = OUTP.matTWISTGLOB(valTIMESTEP-1,:);
    
    [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_t, SURF.matCENTER, COND.valDELTIME, valTIMESTEP);
    
    [SURF, OUTP, INPU] = fcnBEAMFORCE(SURF, OUTP, COND, INPU, FLAG, valTIMESTEP);
    
    if COND.valGUSTTIME == 1
        COND.valDELTIME_old = COND.valDELTIME;
    end
    
    [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME_old,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old,COND.start_loc);



elseif valTIMESTEP >= COND.valSTIFFSTEPS + 1
    
    [COND, INPU, OUTP, MISC, SURF, FLAG, TRIM, VEHI] = fcnFLEXVEHI(INPU, COND, SURF, OUTP, FLAG, MISC, VEHI, TRIM, valTIMESTEP);
    n = n + 1;
    
% elseif FLAG.FLIGHTDYN == 1 && valTIMESTEP >= 30
% 
%     SURF.matCENTER_old = SURF.matCENTER;
%     [SURF, MISC, COND, INPU] = fcnMOVEVEHI(COND, SURF, OUTP, INPU, MISC, FLAG, VEHI, TRIM);
%     [SURF.matUINF] = fcnFLEXUINF(SURF.matCENTER_old, SURF.matCENTER, COND.valDELTIME, 1);
    
else
    
    [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING_STATIC(INPU, VEHI, MISC, COND, SURF, VISC, OUTP, valTIMESTEP);
    
end