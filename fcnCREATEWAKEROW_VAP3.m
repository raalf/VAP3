function [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNTVLST, vecWDVEPANEL, ...
    valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING, vecWDVETRI] = fcnCREATEWAKEROW_VAP3(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, ...
    vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, ...
    matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNTVLST, vecDVEPANEL, ...
    vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVEWING, vecWDVEWING, flagSTEADY, valWSIZE, matNEWWAKEPANEL, vecN, vecWDVETRI)

matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
matNPWAKEGEOM = cat(1, matNPWAKEGEOM, matNPNEWWAKE);
len = length(matNEWWAKE(:,1))*4; 

%% Putting matNEWWAKE into P1/P2/P3/P4 format
% Similar to GenerateDVESTRI
panels = max(vecDVEPANEL);
len = sum(vecN).*2;

vecWDVETRI = [vecWDVETRI; ones(len,1)];

P1          = nan(panels*2,3);
P2          = nan(panels*2,3);
P3          = nan(panels*2,3);
P4          = nan(panels*2,3);

vecEnd      = cumsum(vecN.*repmat(2,length(vecN),1));

for i = 1:panels
    count = 4.*(vecN(i)./2);
    idxStart = vecEnd(i)-count+1;
    idxEnd = vecEnd(i);
    
    panel4corners = reshape(permute(matNEWWAKEPANEL(i,:,:), [2 1 3]), size(matNEWWAKEPANEL(i,:,:), 2), [])';
    [ ~, LE_Left, LE_Right, TE_Left, TE_Right ] = fcnPANEL2DVE(panel4corners, i, vecN./2, ones(size(vecN)));

    TE_Mid = (TE_Right + TE_Left)./2;
    LE_Mid = (LE_Right + LE_Left)./2;
    
    P1(idxStart:4:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count/4,3);
    P2(idxStart:4:idxEnd,:) = reshape(permute(LE_Mid, [2 1 3]),count/4,3);
    P3(idxStart:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);
    P4(idxStart:4:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count/4,3);
    
    P1(idxStart+1:4:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count/4,3);
    P2(idxStart+1:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);
    P3(idxStart+1:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);
    P4(idxStart+1:4:idxEnd,:) = reshape(permute(TE_Left, [2 1 3]),count/4,3);
    
    P1(idxStart+2:4:idxEnd,:) = reshape(permute(LE_Mid, [2 1 3]),count/4,3);
    P2(idxStart+2:4:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count/4,3);
    P3(idxStart+2:4:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count/4,3);
    P4(idxStart+2:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);

    P1(idxStart+3:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);
    P2(idxStart+3:4:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count/4,3);
    P3(idxStart+3:4:idxEnd,:) = reshape(permute(TE_Right, [2 1 3]),count/4,3);
    P4(idxStart+3:4:idxEnd,:) = reshape(permute(TE_Mid, [2 1 3]),count/4,3);
    
    
   % Write DVE CENTER POINT Coordinates
   xsi_vec = (P4(idxStart:idxEnd,:) - P1(idxStart:idxEnd,:) + P3(idxStart:idxEnd,:) - P2(idxStart:idxEnd,:))./2;
   temp_wcenter(idxStart:idxEnd,:) = ((P1(idxStart:idxEnd,:) + P2(idxStart:idxEnd,:))./2) + xsi_vec./2;
    
end

matWCENTER = [matWCENTER; temp_wcenter];


%% Getting wake parameter values from fcnDVECORNER2PARAM
[wdve_eta, vecWDVEHVCRD(end+1:end+len,1), vecWDVEROLL(end+1:end+len,1), vecWDVEPITCH(end+1:end+len,1), vecWDVEYAW(end+1:end+len,1), vecWDVELESWP(end+1:end+len,1), vecWDVEMCSWP(end+1:end+len,1), vecWDVETESWP(end+1:end+len,1), ...
    vecWDVEAREA(end+1:end+len,1), matWDVENORM(end+1:end+len,1:3), ...
    newvertices, newdves, ~] = fcnDVECORNER2PARAM(temp_wcenter, P1, P2, P3, P4);

verticeList = [P1;P2;P3;P4];

if isempty(matWVLST) == 1 % if matWVLST is empty. create the first wake row
    [newvertices,~,newdves] = unique(verticeList,'rows');
    newdves = reshape(newdves,len,4);
    matWDVE(end+1:end+len,1:4) = newdves + length(matWVLST);
    
else % deal with duplicate wake vertices
%     verticeList
    [idx,newdves]=ismember(verticeList,matWVLST,'rows');
%     [idx, newdves] = ismembertol(verticeList, matWVLST, 'ByRows',true);
    [newvertices,~,newrowdves]=unique(verticeList(~idx,:),'rows');
%     [newvertices,~,newrowdves] = uniquetol(verticeList(~idx,:),'ByRows',true);
    
    newdves(~idx) = newrowdves + length(matWVLST(:,1));
    
    
    newdves = reshape(newdves,len,4);
    matWDVE(end+1:end+len,1:4) = newdves;
end

%% Assinging remaining values to wake parameters


matWVLST = cat(1, matWVLST, newvertices);

valWNELE = valWNELE + len;

%% Assigning circulation values to wake DVEs
% K_g = A + ((eta.^2)/3) * C
if flagSTEADY == 1
    vecWKGAM = repmat([reshape(repmat(matCOEFF(vecDVETE>0,1),1,2)',[],1) + ((wdve_eta.^2)./3).*reshape(repmat(matCOEFF(vecDVETE>0,3),1,2)',[],1)], valWNELE/(valWSIZE), 1);
else %unsteady is incorrect
    vecWKGAM(end+1:end+len,1) = [repmat(matCOEFF(vecDVETE>0,1),2,1) + ((wdve_eta.^2)./3).*repmat(matCOEFF(vecDVETE>0,3),2,1)];
end


% matWVLST = cat(1, matWVLST, unique([P1; P2],'rows'));

matWCOEFF = cat(1, matWCOEFF, reshape(repmat(matCOEFF(vecDVETE>0,:),1,2)',3,[],1)');

vecWDVEHVSPN(end+1:end+len,1) = wdve_eta;
vecWDVEPANEL = cat(1, vecWDVEPANEL, reshape(repmat(vecDVEPANEL(vecDVETE>0),1,2)',1,[])');
vecWK = cat(1, vecWK, reshape(repmat(vecK(vecDVETE>0),1,2)',1,[])');

vecWDVEWING = cat(1, vecWDVEWING, reshape(repmat(vecDVEWING(vecDVETE > 0),1,2)',1,[])');

if valWNELE - len == 0
    [ matWADJE, vecWDVESYM, vecWDVETIP, ~, ~ ] = fcnDVEADJT(P1, P2, P3, P4, valWNELE, vecWDVEPANEL, vecSYM );
    valLENWADJE = length(matWADJE(:,1));
else
    
    new_adje_spanwise = [matWADJE(1:valLENWADJE,1) + valWNELE - len matWADJE(1:valLENWADJE,2) matWADJE(1:valLENWADJE,3)+valWNELE-len matWADJE(1:valLENWADJE,4)];
    new_adje_te = [[(valWNELE - len + 2):2:valWNELE]' repmat(3,len/2,1) [(valWNELE - 2*len + 1):2:(valWNELE - len)]' ones(len/2,1)];
    old_adje_le = [new_adje_te(:,3) ones(len/2,1) new_adje_te(:,1) ones(len/2,1)];
    
    % [matWADJE]  DVE# | Local Edge | DVE# | # of Panels This DVE is Touching
    matWADJE = [matWADJE(:,1:4); old_adje_le; new_adje_spanwise; new_adje_te];
    vecWDVESYM = [vecWDVESYM; vecWDVESYM(1:len)];
    vecWDVETIP = [vecWDVETIP; vecWDVETIP(1:len)];
end

end

