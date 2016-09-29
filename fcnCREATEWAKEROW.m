function [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE] ...
    = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, vecWDVEPANEL, vecSYM, valLENWADJE)

matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
matNPWAKEGEOM = cat(1, matNPWAKEGEOM, matNPNEWWAKE);
len = length(matNEWWAKE(:,1));

matWCENTER(end+1:end+len,:) = mean(matNEWWAKE,3);


[vecWDVEHVSPN(end+1:end+len,1), vecWDVEHVCRD(end+1:end+len,1), vecWDVEROLL(end+1:end+len,1), vecWDVEPITCH(end+1:end+len,1), vecWDVEYAW(end+1:end+len,1), vecWDVELESWP(end+1:end+len,1), vecWDVEMCSWP(end+1:end+len,1), vecWDVETESWP(end+1:end+len,1), ...
    vecWDVEAREA(end+1:end+len,1), matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(mean(matNEWWAKE,3), matNEWWAKE(:,:,1), matNEWWAKE(:,:,2), matNEWWAKE(:,:,3), matNEWWAKE(:,:,4));

matWDVE(end+1:end+len,1:4) = newdves + length(matWVLST);

valWNELE = valWNELE + len;

matWVLST = cat(1, matWVLST, newvertices);

matWCOEFF = cat(1, matWCOEFF, matCOEFF(vecDVETE>0,:));

vecWDVEPANEL = cat(1, vecWDVEPANEL, vecDVEPANEL(vecDVETE>0));


if valWNELE - len == 0
    [ matWADJE, vecWDVESYM, vecWDVETIP, vecWDVETE ] = fcnDVEADJT(matNPNEWWAKE(:,:,1), matNPNEWWAKE(:,:,2), matNPNEWWAKE(:,:,3), matNPNEWWAKE(:,:,4), valWNELE, vecWDVEPANEL, vecSYM );
    valLENWADJE = length(matWADJE(:,1));
else
    
    new_adje_spanwise = [matWADJE(1:valLENWADJE,1) + valWNELE - len matWADJE(1:valLENWADJE,2) matWADJE(1:valLENWADJE,3) + valWNELE - len];
    
    new_adje_te = [[(valWNELE - len + 1):valWNELE]' repmat(3,len,1) [(valWNELE - 2*len + 1):(valWNELE - len)]'];
    
    old_adje_le = [new_adje_te(:,3) repmat(1,len,1) new_adje_te(:,1)];
    
    matWADJE = [matWADJE; old_adje_le; new_adje_spanwise; new_adje_te];
    
end

end

