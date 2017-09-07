function [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNTVLST, vecWDVEPANEL, ...
    valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, ...
    vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, ...
    matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNTVLST, vecDVEPANEL, ...
    vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVEWING, vecWDVEWING, flagSTEADY, valWSIZE, vecDVETIP)

matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
matNPWAKEGEOM = cat(1, matNPWAKEGEOM, matNPNEWWAKE);
len = length(matNEWWAKE(:,1));

matWCENTER(end+1:end+len,:) = mean(matNEWWAKE,3);

% Getting wake parameter values from fcnDVECORNER2PARAM
[wdve_eta, vecWDVEHVCRD(end+1:end+len,1), vecWDVEROLL(end+1:end+len,1), vecWDVEPITCH(end+1:end+len,1), vecWDVEYAW(end+1:end+len,1), vecWDVELESWP(end+1:end+len,1), vecWDVEMCSWP(end+1:end+len,1), vecWDVETESWP(end+1:end+len,1), ...
    vecWDVEAREA(end+1:end+len,1), matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(mean(matNEWWAKE,3), matNEWWAKE(:,:,1), matNEWWAKE(:,:,2), matNEWWAKE(:,:,3), matNEWWAKE(:,:,4), vecWDVEWING);

valWNELE = valWNELE + len;

% Assigning circulation values to wake DVEs
% K_g = A + ((eta.^2)/3) * C
if flagSTEADY == 1
    vecWKGAM = repmat([matCOEFF(vecDVETE>0,1) + ((wdve_eta.^2)./3).*matCOEFF(vecDVETE>0,3)], valWNELE/valWSIZE, 1);  
else
    vecWKGAM(end+1:end+len,1) = [matCOEFF(vecDVETE>0,1) + ((wdve_eta.^2)./3).*matCOEFF(vecDVETE>0,3)];
end

% Assinging remaining values to wake parameters
matWDVE(end+1:end+len,1:4) = newdves + length(matWVLST);
matWVLST = cat(1, matWVLST, newvertices);
matWCOEFF = cat(1, matWCOEFF, matCOEFF(vecDVETE>0,:));
vecWDVEHVSPN(end+1:end+len,1) = wdve_eta;
vecWDVEPANEL = cat(1, vecWDVEPANEL, vecDVEPANEL(vecDVETE>0));
vecWK = cat(1, vecWK, vecK(vecDVETE>0));
vecWDVEWING = cat(1, vecWDVEWING, vecDVEWING(vecDVETE > 0));

if valWNELE - len == 0
    [ matWADJE, vecWDVESYM, vecWDVETIP, ~, ~ ] = fcnDVEADJT(matNPNEWWAKE(:,:,1), matNPNEWWAKE(:,:,2), matNPNEWWAKE(:,:,3), matNPNEWWAKE(:,:,4), valWNELE, vecWDVEPANEL, vecSYM );
    valLENWADJE = length(matWADJE(:,1));
else
    
    new_adje_spanwise = [matWADJE(1:valLENWADJE,1) + valWNELE - len matWADJE(1:valLENWADJE,2) matWADJE(1:valLENWADJE,3)+valWNELE-len matWADJE(1:valLENWADJE,4)];
    new_adje_te = [[(valWNELE - len + 1):valWNELE]' repmat(3,len,1) [(valWNELE - 2*len + 1):(valWNELE - len)]' ones(len,1)];
    old_adje_le = [new_adje_te(:,3) ones(len,1) new_adje_te(:,1) ones(len,1)];
    
    % [matWADJE]  DVE# | Local Edge | DVE# | # of Panels This DVE is Touching
    matWADJE = [matWADJE(:,1:4); old_adje_le; new_adje_spanwise; new_adje_te];
    vecWDVESYM = [vecWDVESYM; vecWDVESYM(1:len)];
    vecWDVETIP = [vecWDVETIP; vecWDVETIP(1:len)];
end

end

