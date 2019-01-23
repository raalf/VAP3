function [OUTP] = fcnWINGNFORCE(valTIMESTEP, inddrag, SURF, INPU, COND, OUTP, FLAG, WAKE, VISC)
%% Wing Normal Force
% this routine adds up the DVE's normal forces in order to compute the
% total wing normal forces/density and coefficients based on free stream

% INPUT:

% OUTPUT:
% valCL - Lift coefficient (freestream + induced)
% valCLF - lift coefficient due to freestream
% valCLI - Induced lift coefficient
% valCY - Side force coefficient (freestream + induced)
% valCYF - side force coefficient due to freestream
% valCYI - Induced side force coefficient

valCLI = nan(1,INPU.valVEHICLES,1);
valCL = nan(1,INPU.valVEHICLES,1);
valCLF = nan(1,INPU.valVEHICLES,1);
valCY = nan(1,INPU.valVEHICLES,1);
valCYF = nan(1,INPU.valVEHICLES,1);
valCYI = nan(1,INPU.valVEHICLES,1);
valCDI = nan(1,INPU.valVEHICLES,1);
valE = nan(1,INPU.valVEHICLES,1);

for i = 1:INPU.valVEHICLES
    
    idxvehwing = SURF.vecDVEWING > 0 & SURF.vecDVEVEHICLE == i; %(SURF.vecDVEWING.*SURF.vecDVEVEHICLE == i) > 0;
    
    if any(idxvehwing)

        %     q = 0.5.*sum(abs(vecUINF).^2,2).*valAREA;
        q(i) = 0.5.*mean((sqrt(sum(SURF.matUINF(idxvehwing,:).^2,2)).^2)).*INPU.vecAREA(i);
        
        ntfree = zeros(2,1);
        ntind = zeros(2,1);
        inddragsum = 0;
        
        %sum values from all elements
        ntfree(1) = sum(SURF.vecDVELFREE(idxvehwing));
        ntfree(2) = sum(SURF.vecDVESFREE(idxvehwing));
        
        ntind(1) = sum(SURF.vecDVELIND(idxvehwing));
        ntind(2) = sum(SURF.vecDVESIND(idxvehwing));
        
        inddragsum = sum(inddrag(idxvehwing));
        %double the force if we are using symmetry. This only works with sym for
        %the whole system
        if any(SURF.vecDVESYM) == 1 && ~any(COND.vecVEHBETA(i)) %not sure why beta has to be zero
%             disp('Symmetry is not currently working here! fcnWINGNFORCE');
            ntfree(1) = ntfree(1)*2; %why dont we double the side force?
            ntind(1) = ntind(1)*2;
            ntind(2) = ntind(2)*2;
            inddragsum = inddragsum*2;
        end
        
        %non-dimensionalize
        valCL(i) = (ntfree(1) + ntind(1))/q(i);
        valCLF(i) = ntfree(1)/q(i);
        
        valCY(i) = (ntfree(2) + ntind(2))/q(i);
        valCYF(i) = ntfree(2)/q(i);
        
        valCLI(i) = ntind(1)/q(i);
        valCYI(i) = ntind(2)/q(i);
        
        valCDI(i) = inddragsum/q(i);
        
        AR(i) = (INPU.vecSPAN(i)*INPU.vecSPAN(i))/INPU.vecAREA(i);
        valE(i) = (valCL(i)*valCL(i))/ (pi*AR(i)*valCDI(i));
        
        try
            if FLAG.STIFFWING ~= 1

                [SURF.vecLEDVES,~,~] = find(SURF.vecDVELE(SURF.idxFLEX) > 0);
                
                [OUTP, ~, ~] = fcnVISCOUS(valTIMESTEP, OUTP, COND, VISC, SURF, INPU, FLAG, WAKE, 1);
                [OUTP] = fcnFORCEDIST(SURF, COND, INPU, OUTP, valTIMESTEP);
                
                valDY = 0.5*INPU.vecSPAN/INPU.valNSELE;

                temp_y = (0:valDY:0.5*INPU.vecSPAN)';

                SURF.vecSPANDIST(end) = temp_y(end);
                
                [OUTP.vecLIFTDIST] = interp1(SURF.vecSPANDIST,OUTP.vecLIFTDIST,temp_y);
                [OUTP.vecMOMDIST] = interp1(SURF.vecSPANDIST,OUTP.vecMOMDIST',temp_y);

            else

                OUTP.vecLIFTDIST = [];
                OUTP.vecMOMDIST = [];

            end
        catch
        end
        
    end
    
end

OUTP.vecCL(valTIMESTEP,:) = valCL;
OUTP.vecCLF(valTIMESTEP,:) = valCLF;
OUTP.vecCLI(valTIMESTEP,:) = valCLI;
OUTP.vecCY(valTIMESTEP,:) = valCY;
OUTP.vecCYF(valTIMESTEP,:) = valCYF;
OUTP.vecCYI(valTIMESTEP,:) = valCYI;
OUTP.vecCDI(valTIMESTEP,:) = valCDI;
OUTP.vecE(valTIMESTEP,:) = valE;

end