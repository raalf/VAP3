function [vecCLv, vecCD, vecPREQ, vecLD, valCL, valCLF, valCLI, valCDI, valCT, valCP, valE, vecDVENFREE, vecDVENIND, vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecDVESYM, vecDVETESWP, valAREA, valSPAN, valBETA, vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEVEHICLE, valVEHICLES, matVEHROT, flagTRI, flagSTEADY, flagGPU, vecDVEROTOR, matROTORAXIS, vecROTORRPM, vecROTDIAM, valDELTIME, vecVEHVINF, valDENSITY, valKINV,  vecDVEAREA, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, valFPWIDTH, valINTERF, valMAXTIME, vecPANELWING, matUINFROT, vecQARM, flagVISCOUS)
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

% %% Viscous wrapper
[vecCLv, vecCD, vecPREQ, vecLD, matROTORCDP] = fcnVISCOUS(valCL, valCDI, vecVEHVINF, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
    matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
    valFPWIDTH, valINTERF, vecDVEROLL, valVEHICLES, vecDVEVEHICLE, vecDVEROTOR, matUINF, valTIMESTEP, valMAXTIME, valDELTIME, vecROTORRPM, matDVE, matVLST,...
    matCOEFF, vecK, vecDVEHVSPN, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, flagTRI, flagSTEADY, flagGPU, valNELE, vecPANELWING, flagVISCOUS);

%% Rotor Forces
if max(vecDVEROTOR) > 0
    [valCT, valCP] = fcnROTORFORCE(matROTORCDP, en, vecDVENFREE, vecDVENIND, inddrag, matUINF, vecDVEROTOR, matVEHROT, matROTORAXIS, vecROTORRPM, vecROTDIAM, matUINFROT, vecQARM);
else
    valCT = nan;
    valCP = nan;
end

end

