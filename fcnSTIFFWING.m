function [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valTIMESTEP)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

[SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);

n = 1;

% if valGUSTTIME > 1 || valTIMESTEP == valGUSTSTART
    
%     [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old);
%     valGUSTTIME = valGUSTTIME + 1;
    
% end

OUTP.matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN(FLAG.vecFLEXIBLE == 1),1)+1);

OUTP.matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN(FLAG.vecFLEXIBLE == 1),1)+1);

end