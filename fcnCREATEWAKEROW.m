function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF] = fcnCREATEWAKEROW(FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF)

WAKE.matWAKEGEOM = cat(1, WAKE.matWAKEGEOM, MISC.matNEWWAKE);
WAKE.matNPWAKEGEOM = cat(1, WAKE.matNPWAKEGEOM, MISC.matNPNEWWAKE);
len = length(MISC.matNEWWAKE(:,1));

WAKE.matWCENTER(end+1:end+len,:) = mean(MISC.matNEWWAKE,3);

% Getting wake parameter values from fcnDVECORNER2PARAM
[wdve_eta, WAKE.vecWDVEHVCRD(end+1:end+len,1), WAKE.vecWDVEROLL(end+1:end+len,1), WAKE.vecWDVEPITCH(end+1:end+len,1), WAKE.vecWDVEYAW(end+1:end+len,1), WAKE.vecWDVELESWP(end+1:end+len,1), WAKE.vecWDVEMCSWP(end+1:end+len,1), WAKE.vecWDVETESWP(end+1:end+len,1), ...
    WAKE.vecWDVEAREA(end+1:end+len,1), WAKE.matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(mean(MISC.matNEWWAKE,3), MISC.matNEWWAKE(:,:,1), MISC.matNEWWAKE(:,:,2), MISC.matNEWWAKE(:,:,3), MISC.matNEWWAKE(:,:,4), WAKE.vecWDVESURFACE);

WAKE.valWNELE = WAKE.valWNELE + len;

% Assigning circulation values to wake DVEs
% K_g = A + ((eta.^2)/3) * C
if FLAG.STEADY == 1
    WAKE.vecWKGAM = repmat([SURF.matCOEFF(SURF.vecDVETE>0,1) + ((wdve_eta.^2)./3).*SURF.matCOEFF(SURF.vecDVETE>0,3)], WAKE.valWNELE/WAKE.valWSIZE, 1);  
else
    WAKE.vecWKGAM(end+1:end+len,1) = [SURF.matCOEFF(SURF.vecDVETE>0,1) + ((wdve_eta.^2)./3).*SURF.matCOEFF(SURF.vecDVETE>0,3)];
end

% Assinging remaining values to wake parameters
WAKE.matWDVE(end+1:end+len,1:4) = newdves + length(WAKE.matWVLST);
WAKE.matWVLST = cat(1, WAKE.matWVLST, newvertices);
WAKE.matWCOEFF = cat(1, WAKE.matWCOEFF, SURF.matCOEFF(SURF.vecDVETE>0,:));
WAKE.vecWDVEHVSPN(end+1:end+len,1) = wdve_eta;
WAKE.vecWDVEPANEL = cat(1, WAKE.vecWDVEPANEL, SURF.vecDVEPANEL(SURF.vecDVETE>0));
WAKE.vecWK = cat(1, WAKE.vecWK, SURF.vecK(SURF.vecDVETE>0));
WAKE.vecWDVESURFACE = cat(1, WAKE.vecWDVESURFACE, SURF.vecDVESURFACE(SURF.vecDVETE > 0));
WAKE.vecWPLOTSURF = cat(1, WAKE.vecWPLOTSURF, SURF.vecDVEWING(SURF.vecDVETE > 0) + (SURF.vecDVEROTOR(SURF.vecDVETE > 0) + max(SURF.vecDVEWING).*uint8(SURF.vecDVEROTOR(SURF.vecDVETE > 0) > 0)));

if WAKE.valWNELE - len == 0
    [ WAKE.matWADJE, WAKE.vecWDVESYM, WAKE.vecWDVETIP, ~, ~ ] = fcnDVEADJT(MISC.matNPNEWWAKE(:,:,1), MISC.matNPNEWWAKE(:,:,2), MISC.matNPNEWWAKE(:,:,3), MISC.matNPNEWWAKE(:,:,4), WAKE.valWNELE, WAKE.vecWDVEPANEL, INPU.vecSYM );
    WAKE.valLENWADJE = length(WAKE.matWADJE(:,1));
else
    
    new_adje_spanwise = [WAKE.matWADJE(1:WAKE.valLENWADJE,1) + WAKE.valWNELE - len WAKE.matWADJE(1:WAKE.valLENWADJE,2) WAKE.matWADJE(1:WAKE.valLENWADJE,3)+WAKE.valWNELE-len WAKE.matWADJE(1:WAKE.valLENWADJE,4)];
    new_adje_te = [[(WAKE.valWNELE - len + 1):WAKE.valWNELE]' repmat(3,len,1) [(WAKE.valWNELE - 2*len + 1):(WAKE.valWNELE - len)]' ones(len,1)];
    old_adje_le = [new_adje_te(:,3) ones(len,1) new_adje_te(:,1) ones(len,1)];
    
    % [WAKE.matWADJE]  DVE# | Local Edge | DVE# | # of Panels This DVE is Touching
    WAKE.matWADJE = uint32([WAKE.matWADJE(:,1:4); old_adje_le; new_adje_spanwise; new_adje_te]);
    WAKE.vecWDVESYM = [WAKE.vecWDVESYM; WAKE.vecWDVESYM(1:len)];
    WAKE.vecWDVETIP = [WAKE.vecWDVETIP; WAKE.vecWDVETIP(1:len)];
end

end

