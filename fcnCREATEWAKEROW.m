function [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM] ...
    = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK)

matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
matNPWAKEGEOM = cat(1, matNPWAKEGEOM, matNPNEWWAKE);
len = length(matNEWWAKE(:,1));

matWCENTER(end+1:end+len,:) = mean(matNEWWAKE,3);

% Getting wake parameter values from fcnDVECORNER2PARAM
[vecWDVEHVSPN(end+1:end+len,1), wdve_eta, vecWDVEROLL(end+1:end+len,1), vecWDVEPITCH(end+1:end+len,1), vecWDVEYAW(end+1:end+len,1), vecWDVELESWP(end+1:end+len,1), vecWDVEMCSWP(end+1:end+len,1), vecWDVETESWP(end+1:end+len,1), ...
    vecWDVEAREA(end+1:end+len,1), matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(mean(matNEWWAKE,3), matNEWWAKE(:,:,1), matNEWWAKE(:,:,2), matNEWWAKE(:,:,3), matNEWWAKE(:,:,4));

% Assigning circulation values to wake DVEs
% K_g = A + ((eta.^2)/3) * C
vecWKGAM(end+1:end+len,1) = [matCOEFF(vecDVETE>0,1) + ((wdve_eta.^2)./3).*matCOEFF(vecDVETE>0,3)];

% Assinging remaining values to wake parameters
matWDVE(end+1:end+len,1:4) = newdves + length(matWVLST);
valWNELE = valWNELE + len;
matWVLST = cat(1, matWVLST, newvertices);
matWCOEFF = cat(1, matWCOEFF, matCOEFF(vecDVETE>0,:));
vecWDVEHVCRD(end+1:end+len,1) = wdve_eta;
vecWDVEPANEL = cat(1, vecWDVEPANEL, vecDVEPANEL(vecDVETE>0));
vecWK = cat(1, vecWK, vecK);

if valWNELE - len == 0
    [ matWADJE, vecWDVESYM, vecWDVETIP, ~, ~ ] = fcnDVEADJT(matNPNEWWAKE(:,:,1), matNPNEWWAKE(:,:,2), matNPNEWWAKE(:,:,3), matNPNEWWAKE(:,:,4), valWNELE, vecWDVEPANEL, vecSYM );
    valLENWADJE = length(matWADJE(:,1));
else
    
    new_adje_spanwise = [matWADJE(1:valLENWADJE,1) + valWNELE - len matWADJE(1:valLENWADJE,2) matWADJE(1:valLENWADJE,3) + valWNELE - len];
    
    new_adje_te = [[(valWNELE - len + 1):valWNELE]' repmat(3,len,1) [(valWNELE - 2*len + 1):(valWNELE - len)]'];
    
    old_adje_le = [new_adje_te(:,3) repmat(1,len,1) new_adje_te(:,1)];
    
    matWADJE = [matWADJE; old_adje_le; new_adje_spanwise; new_adje_te];
    vecWDVESYM = [vecWDVESYM; vecWDVESYM(1:len)];
    vecWDVETIP = [vecWDVETIP; vecWDVETIP(1:len)];
end

end

