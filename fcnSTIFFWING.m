function [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valTIMESTEP)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

[SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);

if valGUSTTIME > 1 || valTIMESTEP == valGUSTSTART
        [matUINF, gust_vel] = fcnGUSTWINGSTIFF(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valUINF,valGUSTSTART,fpg,gust_vel_old)
    valGUSTTIME = valGUSTTIME + 1;
end
OUTP.matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN,1)+1);
OUTP.matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN,1)+1);

end