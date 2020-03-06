function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMVEH(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, ii)

% Find alpha derivatives
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnALPHASWEEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

% Find elevator derivatives
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnELEVSWEEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

trim = [TRIM.CLalpha*pi/180, TRIM.CLde*pi/180; TRIM.Cmalpha*pi/180, TRIM.Cmde*pi/180]\[COND.CLtrim; -TRIM.Cm0];

% Compute how much to rotate alpha and tail to be in trimmed state
trim_rot(1) = trim(1)+TRIM.alpha0 - COND.vecVEHALPHA;
trim_rot(2) = deg2rad(TRIM.tau*(trim(2)-SURF.vecELEVANGLE*180/pi));

% Store alpha and elevator angle for trim at each iteration for comparison
% later
OUTP.valALPHATRIM(ii,1) = trim(1)+TRIM.alpha0;
OUTP.valTAILTRIM(ii,1) = trim(2);

SURF.vecELEVANGLE = OUTP.valTAILTRIM(ii);

COND.vecVEHALPHA = OUTP.valALPHATRIM(ii);

[INPU, VEHI, SURF, MISC] = fcnSETTRIM(FLAG, COND, INPU, VEHI, SURF, MISC, trim_rot);

end