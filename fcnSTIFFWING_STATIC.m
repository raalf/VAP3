function [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING_STATIC(INPU, VEHI, MISC, COND, SURF, VISC, OUTP, valTIMESTEP)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.


[SURF, INPU, MISC, VISC] = fcnMOVESURFACE_STATIC(INPU, VEHI, MISC, COND, SURF, VISC);

if valTIMESTEP == 1
    OUTP.matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN,1)+1);
    OUTP.matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(INPU.vecN,1)+1);
else
    OUTP.matDEFGLOB(valTIMESTEP,:) = OUTP.matDEFGLOB(valTIMESTEP - 1,:);
    OUTP.matTWISTGLOB(valTIMESTEP,:) = OUTP.matTWISTGLOB(valTIMESTEP - 1,:);
end

end