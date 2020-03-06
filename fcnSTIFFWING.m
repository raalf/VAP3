function [SURF, INPU, MISC, VISC, OUTP, VEHI] = fcnSTIFFWING(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valTIMESTEP)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

[SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);

if any(FLAG.vecTRIMABLE == 1) == 1
    
    SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) + VEHI.matVEHUVW*COND.valDELTIME;
    
end

if any(FLAG.vecFLEXIBLE == 1) == 1
   
    SURF.matEALST = SURF.matEALST + VEHI.matVEHUVW*COND.valDELTIME;
    SURF.vecWINGCG = SURF.vecWINGCG + VEHI.matVEHUVW*COND.valDELTIME;
    VEHI.vecFUSECG = VEHI.vecFUSECG + VEHI.matVEHUVW*COND.valDELTIME;
    VEHI.vecWINGCG(2:end,:) = VEHI.vecWINGCG(2:end,:) + VEHI.matVEHUVW*COND.valDELTIME;
    VEHI.vecPAYLCG = VEHI.vecPAYLCG + VEHI.matVEHUVW*COND.valDELTIME;
    VEHI.vecPROPLOC = VEHI.vecPROPLOC + VEHI.matVEHUVW*COND.valDELTIME;
    
end

n = 1;

% if valGUSTTIME > 1 || valTIMESTEP == valGUSTSTART
    
%     [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old);
%     valGUSTTIME = valGUSTTIME + 1;
    
% end

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

end