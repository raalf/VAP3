function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient
% valCDI - Ind. drag coefficient
% valE - Span Efficiency

%% Element normal forces, lift forces and side forces (freestream and induced)
[en, SURF.vecDVENFREE, SURF.vecDVENIND, SURF.vecDVELFREE, SURF.vecDVELIND, SURF.vecDVESFREE, SURF.vecDVESIND] = fcnDVENFORCE(valTIMESTEP, SURF, WAKE, VEHI, FLAG, INPU);

%% Induced Drag force
[inddrag] = fcnDVEINDDRAG(valTIMESTEP, SURF, WAKE, FLAG);

%% Sum up element forces to generate total wing forces
OUTP = fcnWINGNFORCE(valTIMESTEP, inddrag, SURF, INPU, COND, OUTP);

%% Viscous wrapper
[OUTP, matROTORCDP, vecDELNDIST] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE);

%% Rotor Forces
if max(SURF.vecDVEROTOR) > 0
    [OUTP.vecCT(valTIMESTEP,:), OUTP.vecCP(valTIMESTEP,:)] = fcnROTORFORCE(matROTORCDP, en, SURF.vecDVENFREE, SURF.vecDVENIND, inddrag, SURF.matUINF, SURF.vecDVEROTOR, VEHI.matVEHROT, INPU.matROTORAXIS, COND.vecROTORRPM, INPU.vecROTDIAM, SURF.matUINFROT, SURF.vecQARM, SURF.vecDVEVEHICLE, vecDELNDIST);
else
    OUTP.vecCT = nan;
    OUTP.vecCP = nan;
end

end

