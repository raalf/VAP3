function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP, FLAG] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient
% valCDI - Ind. drag coefficient
% valE - Span Efficiency

%% Element normal forces, lift forces and side forces (freestream and induced)
[en, SURF.vecDVENFREE, SURF.vecDVENIND, SURF.vecDVELFREE, SURF.vecDVELIND, SURF.vecDVESFREE, SURF.vecDVESIND, SURF.matLIFTDIR] = fcnDVENFORCE(valTIMESTEP, SURF, WAKE, VEHI, FLAG, INPU);

%% Induced Drag force
[inddrag] = fcnDVEINDDRAG(valTIMESTEP, SURF, WAKE, FLAG);
% OUTP.WING.vecDIDIST(valTIMESTEP,:) = inddrag(SURF.vecDVEWING > 0)';
% ^ JCole shit. Commented out for now. Fix in San Diego

%% Sum up element forces to generate total wing forces
OUTP = fcnWINGNFORCE(valTIMESTEP, inddrag, SURF, INPU, COND, OUTP, FLAG, WAKE, VISC);

%% Viscous wrapper
[OUTP, matROTORCDP, vecDELNDIST] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE, 0);

%% Perform longitudinal trim routine
if valTIMESTEP == COND.valMAXTIME && any(FLAG.vecTRIMABLE == 1) == 1 && FLAG.TRIM == 1
    [OUTP, SURF, FLAG] = fcnLONGTRIM(SURF, INPU, OUTP, FLAG, COND);
end

%% Rotor Forces
if max(SURF.vecDVEROTOR) > 0
    OUTP = fcnROTORFORCE(valTIMESTEP, matROTORCDP, en, inddrag, vecDELNDIST, SURF, VEHI, INPU, COND, OUTP);
else
    OUTP.vecCT = nan;
    OUTP.vecCP = nan;
    OUTP.vecCPI = nan;
end

end

