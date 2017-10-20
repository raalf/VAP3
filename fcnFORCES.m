function [vecCTCONV, valCL, valCLF, valCLI, valCDI, valCT, valE, vecDVENFREE, vecDVENIND, vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecDVESYM, vecDVETESWP, valAREA, valSPAN, valBETA, vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEVEHICLE, valVEHICLES, matVEHROT, flagTRI, flagSTEADY, flagGPU, vecDVEROTOR, matROTORAXIS, vecROTORRPM, vecROTDIAM, valDELTIME, vecCTCONV)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient
% valCDI - Ind. drag coefficient
% valE - Span Efficiency

%% Element normal forces, lift forces and side forces (freestream and induced)

[en, vecDVENFREE, vecDVENIND, vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnDVENFORCE(matCOEFF...
    ,vecK,matDVE,valNELE,matCENTER,matVLST,matUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEHVCRD,vecDVEROLL,...
    vecDVEPITCH,vecDVEYAW,vecDVELE,matADJE,valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD,...
    vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecDVESYM, vecDVETESWP, matVEHROT, vecDVEVEHICLE, flagTRI, flagSTEADY, flagGPU);
%% TE Element induced drag forces

%% Induced Drag force

[inddrag] = fcnDVEINDDRAG(valNELE, matCOEFF, matDVE, matVLST, matUINF, vecDVEHVSPN, vecDVEHVCRD, vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, ...
    valWSIZE, valTIMESTEP, vecDVESYM, vecDVEWING, vecWDVEWING, flagTRI, flagSTEADY, flagGPU);

%% Sum up element forces to generate total wing forces

[valCL, valCLF, valCLI, valCY, valCYF, valCYI, valCDI, valE]= fcnWINGNFORCE(vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, inddrag, matUINF, valAREA, valSPAN, vecDVESYM, valBETA, vecDVEVEHICLE, vecDVEWING, valVEHICLES);

if max(vecDVEROTOR) > 0
    [valCT, vecCTCONV] = fcnROTORFORCE( en, vecDVENFREE, vecDVENIND, inddrag, matUINF, vecDVEROTOR, matVEHROT, matROTORAXIS, vecROTORRPM, vecROTDIAM, matCENTER, valDELTIME, valTIMESTEP, vecCTCONV);
else
    valCT = nan;
end

end

