function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP)
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

%% Sum up element forces to generate total wing forces
OUTP = fcnWINGNFORCE(valTIMESTEP, inddrag, SURF, INPU, COND, OUTP, FLAG, WAKE, VISC);

if ~any(SURF.vecDVEROTOR)
    [OUTP, SURF] = fcnSURFACEDIST(COND, INPU, SURF, OUTP, valTIMESTEP);
end

%% Viscous wrapper
[OUTP, matROTORCDP, vecDELNDIST] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE, 0);

%% Rotor Forces
if max(SURF.vecDVEROTOR) > 0
    OUTP = fcnROTORFORCE(valTIMESTEP, matROTORCDP, en, inddrag, vecDELNDIST, SURF, VEHI, INPU, COND, OUTP);
else
    OUTP.vecCT = nan;
    OUTP.vecCP = nan;
    OUTP.vecCPI = nan;
end

end

