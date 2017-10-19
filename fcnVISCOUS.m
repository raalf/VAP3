function [vecCLv, vecCD, vecPREQ, vecVINF, vecLD] = fcnVISCOUS(vecCL, vecCDI, valWEIGHT, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
    matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
    valFPWIDTH, valINTERF, vecDVEROLL, valVEHICLES, vecDVEVEHICLE)

for i = 1:valVEHICLES
    
    idxvehwing = vecDVEWING > 0 & vecDVEVEHICLE == i; %(vecDVEWING.*vecDVEVEHICLE == i) > 0;
    
    if any(idxvehwing)

        [vecCLv(i), vecCD(i), vecPREQ(i), vecVINF(i), vecLD(i)] = fcnVISCOUS_WING(vecCL(end), vecCDI(end), ...
            valWEIGHT, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
            vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING.*idxvehwing, vecN, vecM, vecDVEAREA, ...
            matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
            valFPWIDTH, valINTERF, vecDVEROLL);

    end
    
    
    idxvehrotor = vecDVEROTOR > 0 & vecDVEVEHICLE == i;
    
    if any(idxvehrotor)
        
       
        
    end
    
end  

end

