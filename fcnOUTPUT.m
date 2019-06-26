function [OUTP] = fcnOUTPUT(COND, FLAG, SURF, OUTP, valTIMESTEP)

% Setting up variables for outputting. Includes time-averaged forces (if
% applicable).

if FLAG.PREVIEW ~= 1 && max(SURF.vecDVEROTOR) > 0 && ~isempty(valTIMESTEP)
    % Time averaged lift
    OUTP.vecCL_AVG = fcnTIMEAVERAGE(OUTP.vecCLv, COND.vecROTORRPM, COND.valDELTIME);
    
    % Time averaged drags (total, induced, profile)
    OUTP.vecCD_AVG = fcnTIMEAVERAGE(OUTP.vecCD, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCDI_AVG = fcnTIMEAVERAGE(OUTP.vecCDI, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCDP_AVG = fcnTIMEAVERAGE(OUTP.vecCD - OUTP.vecCDI, COND.vecROTORRPM, COND.valDELTIME);
   
    OUTP.vecCT_AVG = fcnTIMEAVERAGE(OUTP.vecCT, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCFx_AVG = fcnTIMEAVERAGE(OUTP.vecCFx, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCFy_AVG = fcnTIMEAVERAGE(OUTP.vecCFy, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCMx_AVG = fcnTIMEAVERAGE(OUTP.vecCMx, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCMy_AVG = fcnTIMEAVERAGE(OUTP.vecCMy, COND.vecROTORRPM, COND.valDELTIME);
    
    % Time averaged propeller powers
    OUTP.vecCP_AVG = fcnTIMEAVERAGE(OUTP.vecCP, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCPI_AVG = fcnTIMEAVERAGE(OUTP.vecCPI, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCPP_AVG = fcnTIMEAVERAGE(OUTP.vecCP - OUTP.vecCPI, COND.vecROTORRPM, COND.valDELTIME);

    for i = 1:max(SURF.vecDVEROTOR)
        OUTP.ROTOR(i).vecTHRUSTDIST(any(~any(OUTP.ROTOR(i).vecTHRUSTDIST, 2), 3), :) = [];
        OUTP.ROTOR(i).vecTHRUSTDIST_AVG = fcnTIMEAVERAGE(OUTP.ROTOR(i).vecTHRUSTDIST, COND.vecROTORRPM, COND.valDELTIME);
        OUTP.ROTOR(i).vecTORQUEDIST(any(~any(OUTP.ROTOR(i).vecTORQUEDIST, 2), 3), :) = [];
        OUTP.ROTOR(i).vecTORQUEDIST_AVG = fcnTIMEAVERAGE(OUTP.ROTOR(i).vecTORQUEDIST, COND.vecROTORRPM, COND.valDELTIME);
    end
    
elseif FLAG.PREVIEW ~= 1 && ~any(SURF.vecDVEROTOR)
    OUTP.vecCLv_AVG = OUTP.vecCLv(end);
    OUTP.vecCD_AVG = OUTP.vecCD(end);
    OUTP.vecCDI_AVG = OUTP.vecCDI(end);
    OUTP.vecCDP_AVG = OUTP.vecCD(end) - OUTP.vecCDI(end);
end

OUTP.vecVEHALPHA = COND.vecVEHALPHA;
OUTP.vecCOLLECTIVE = COND.vecCOLLECTIVE;
OUTP.vecROTDIAM = INPU.vecROTDIAM;
OUTP.vecVEHWEIGHT = COND.vecVEHWEIGHT;
OUTP.vecROTORRPM = COND.vecROTORRPM;
OUTP.vecDVEAREA = SURF.vecDVEAREA;
OUTP.valAREA = INPU.vecAREA;
OUTP.matGEOM = INPU.matGEOM;
OUTP.valDENSITY = COND.valDENSITY;