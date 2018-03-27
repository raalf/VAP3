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
%% TE Element induced drag forces

%% Induced Drag force

[inddrag] = fcnDVEINDDRAG(SURF.valNELE, SURF.matCOEFF, SURF.matDVE, SURF.matVLST, SURF.matUINF, SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVETE, WAKE.valWNELE, WAKE.matWDVE, WAKE.matWVLST, WAKE.matWCOEFF, WAKE.vecWK, WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD, WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, WAKE.vecWDVETESWP, ...
    WAKE.valWSIZE, valTIMESTEP, SURF.vecDVESYM, SURF.vecDVEWING, WAKE.vecWDVESURFACE, FLAG.TRI, FLAG.STEADY, FLAG.GPU);

%% Sum up element forces to generate total wing forces
[OUTP.vecCL(valTIMESTEP,:), OUTP.vecCLF(valTIMESTEP,:), OUTP.vecCLI(valTIMESTEP,:), OUTP.vecCY(valTIMESTEP,:), OUTP.vecCYF(valTIMESTEP,:), OUTP.vecCYI(valTIMESTEP,:), OUTP.vecCDI(valTIMESTEP,:), OUTP.vecE(valTIMESTEP,:)]= fcnWINGNFORCE(SURF.vecDVELFREE, SURF.vecDVELIND, SURF.vecDVESFREE, SURF.vecDVESIND, inddrag, SURF.matUINF, INPU.vecAREA, INPU.vecSPAN, SURF.vecDVESYM, COND.vecVEHBETA, SURF.vecDVEVEHICLE, SURF.vecDVEWING, INPU.valVEHICLES);

%% Viscous wrapper
[OUTP, matROTORCDP, vecDELNDIST] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE);

%% Rotor Forces
if max(SURF.vecDVEROTOR) > 0
    [OUTP.vecCT(valTIMESTEP,:), OUTP.vecCP(valTIMESTEP,:)] = fcnROTORFORCE(matROTORCDP, en, SURF.vecDVENFREE, SURF.vecDVENIND, inddrag, SURF.matUINF, SURF.vecDVEROTOR, VEHI.matVEHROT, INPU.matROTORAXIS, COND.vecROTORRPM, INPU.vecROTDIAM, SURF.matUINFROT, SURF.vecQARM, SURF.vecDVEVEHICLE, vecDELNDIST);
else
    OUTP.vecCT = nan;
    OUTP.vecCP = nan;
end

OUTP.vecCTCONV = OUTP.vecCT;

end

