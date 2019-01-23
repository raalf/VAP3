function [OUTP, SURF, FLAG] = fcnLONGTRIM(SURF, INPU, OUTP, FLAG, COND)
% Function to compute new tail incidence angle to trim the aircraft

[SURF, OUTP] = fcnPITCHMOMENT(FLAG, SURF, OUTP, INPU, COND);

for i = 1:INPU.valVEHICLES
    
    if COND.valTRIMSTEP == 1
        if OUTP.vecVEHCM(i) > 0
            deltaEPS = 2*pi/180; % Deflect tail TE 2 deg down
            SURF.vecTAILPITCHold = SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180;
            OUTP.vecVEHCMold(i) = OUTP.vecVEHCM(i);
            
            [SURF] = fcnTAILTRIM(SURF, FLAG, COND, deltaEPS, i);
            
        else
            deltaEPS = -2*pi/180; % Deflect tail TE 2 deg up
            SURF.vecTAILPITCHold = SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180;
            OUTP.vecVEHCMold(i) = OUTP.vecVEHCM(i);
            
            [SURF] = fcnTAILTRIM(SURF, FLAG, COND, deltaEPS, i);
        end
        
    else
        
        tol = OUTP.vecVEHCM(i)*OUTP.vecVEHCM(i);

        if tol < 0.00000025
            FLAG.TRIMMED(i) = 1;
            break;
        end
    
        % Compute new tail incidence angle per vehicle
        tempS = (SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180) -...
            OUTP.vecVEHCM*((SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180)...
            - SURF.vecTAILPITCHold)/(OUTP.vecVEHCM(i)-OUTP.vecVEHCMold(i));
        
        OUTP.vecVEHCMold(i) = OUTP.vecVEHCM(i);
        
        deltaEPS = tempS - (SURF.vecDVEPITCH(SURF.vecDVEWING == 2)-COND.vecVEHALPHA(i)*pi/180);
        
        [SURF] = fcnTAILTRIM(SURF, FLAG, COND, deltaEPS(1), i);
    end
    
end

end