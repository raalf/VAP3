function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP, FLAG] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient
% valCDI - Ind. drag coefficient
% valE - Span Efficiency

%% Element normal forces, lift forces and side forces (freestream and induced)
[en, SURF.vecDVENFREE, SURF.vecDVENIND, SURF.vecDVELFREE, SURF.vecDVELIND, SURF.vecDVESFREE, SURF.vecDVESIND, SURF.matLIFTDIR, SURF.gamma_old, SURF.dGammadt, SURF.wake_vel_time] = fcnDVENFORCE(valTIMESTEP, COND, SURF, WAKE, VEHI, FLAG, INPU);

SURF.nfree(:,valTIMESTEP) = SURF.vecDVENFREE;
SURF.gammaold(:,valTIMESTEP) = SURF.gamma_old;
SURF.en_t(:,:,valTIMESTEP) = en;

%% Induced Drag force
[inddrag, tempUINF] = fcnDVEINDDRAG(valTIMESTEP, SURF, WAKE, FLAG);

SURF.vecDVEDIND = inddrag;
SURF.matDRAGDIR = tempUINF;

%% Sum up element forces to generate total wing forces
[OUTP, SURF] = fcnWINGNFORCE(valTIMESTEP, inddrag, SURF, INPU, COND, OUTP, FLAG, WAKE, VISC);

if ~any(SURF.vecDVEROTOR)
    [OUTP, SURF] = fcnSURFACEDIST(COND, INPU, SURF, OUTP, FLAG, valTIMESTEP);
end

%% Viscous wrapper
[OUTP, SURF, matROTORCDP, vecDELNDIST] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE, 0);

%% Moment calculations
% If there's an h-stab, compute pitch moment
if any(SURF.vecWINGTYPE == 2)
[SURF, OUTP] = fcnPITCHMOMENT(FLAG, SURF, OUTP, INPU, COND, VEHI);
end

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

