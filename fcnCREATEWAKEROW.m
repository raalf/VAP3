function   [matWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER] ...
    = fcnCREATEWAKEROW(matNEWWAKE, matWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER)

matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
len = length(matNEWWAKE(:,1));

matWCENTER(end+1:end+len,:) = mean(matNEWWAKE,3);


[vecWDVEHVSPN(end+1:end+len,1), vecWDVEHVCRD(end+1:end+len,1), vecWDVEROLL(end+1:end+len,1), vecWDVEPITCH(end+1:end+len,1), vecWDVEYAW(end+1:end+len,1), vecWDVELESWP(end+1:end+len,1), vecWDVEMCSWP(end+1:end+len,1), vecWDVETESWP(end+1:end+len,1), ...
    vecWDVEAREA(end+1:end+len,1), matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(mean(matNEWWAKE,3), matNEWWAKE(:,:,1), matNEWWAKE(:,:,2), matNEWWAKE(:,:,3), matNEWWAKE(:,:,4));

matWDVE(end+1:end+len,1:4) = newdves + length(matWVLST);

valWNELE = valWNELE + len;

matWVLST = cat(1, matWVLST, newvertices);

end

