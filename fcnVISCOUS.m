function [vecCLv, vecCD, vecPREQ, vecLD] = fcnVISCOUS(vecCL, vecCDI, vecVEHVINF, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
    matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
    valFPWIDTH, valINTERF, vecDVEROLL, valVEHICLES, vecDVEVEHICLE, vecDVEROTOR, matUINF)

vecCLv = nan(valVEHICLES,1);
vecCD = nan(valVEHICLES,1);
vecPREQ = nan(valVEHICLES,1);
vecLD = nan(valVEHICLES,1);

for i = 1:valVEHICLES
    
    idxvehwing = vecDVEWING > 0 & vecDVEVEHICLE == i; %(vecDVEWING.*vecDVEVEHICLE == i) > 0;
    
    if any(idxvehwing)

        [vecCLv(i), vecCD(i), vecPREQ(i), vecLD(i)] = fcnVISCOUS_WING(vecCL(end), vecCDI(end), ...
            vecVEHVINF(i), valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
            vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING.*idxvehwing, vecN, vecM, vecDVEAREA, ...
            matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
            valFPWIDTH, valINTERF, vecDVEROLL, matUINF);

    end
    
    
    idxvehrotor = vecDVEROTOR > 0 & vecDVEVEHICLE == i;
    
    if any(idxvehrotor)
        
       
        
    end
    
end  

end

