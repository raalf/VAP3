function [ matVLST0, matCENTER0, matDVE, matADJE, vecDVEVEHICLE, ...
    vecDVEWING, vecDVEROTOR, matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, ...
    vecDVETIP, vecDVELE, vecDVETE, vecDVEROTORBLADE, vecDVESYM, ...
    valNELE, matNTVLST0] = fcnDUPBLADE( vecROTORVEH, vecDVEROTOR, ...
    matVLST0, matCENTER0, matDVE, matADJE, vecROTORBLADES, ...
    valNELE, matROTORHUB, matVEHORIG, vecDVEVEHICLE, vecDVEWING, ...
    matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, vecDVETIP, vecDVELE, ...
    vecDVETE, vecDVEROTORBLADE, vecDVESYM, matROTORAXIS, matNTVLST0)

%FCNDUPBLADE Summary of this function goes here
%   Duplicate blades within a rotor
% Add Duplicate Rotor Blades



valROTORS = length(vecROTORVEH);
for n = 1:valROTORS
    idxDVEBLADE = find(vecDVEROTOR==n);
    matDVEBLADE = matDVE(idxDVEBLADE,:);
    [idxVLSTBLADE,~,c] = unique(matDVEBLADE);
    offsetDVEBLADE = reshape(c,[],4);
    
    % matADJE
    checkADJE = ismember(matADJE,idxDVEBLADE);
    offsetADJE = matADJE(and(checkADJE(:,1),checkADJE(:,3)),:);
    
    for j = 1:vecROTORBLADES(n)-1
        
        radBLADE = 2*pi/vecROTORBLADES(n)*j; % radian
        dcmBLADE = angle2dcm(radBLADE,0,0, 'ZXY');
        
        %matDVE
        addDVEBLADE = offsetDVEBLADE + length(matVLST0(:,1));
        matDVE = [matDVE;addDVEBLADE];
        
        %matADJE
        offsetADJEBLADE = offsetADJE;
        offsetADJEBLADE(:,[1,3]) = offsetADJE(:,[1,3]) - min(idxDVEBLADE) + 1 + valNELE;
        matADJE = [matADJE;offsetADJEBLADE];
        
        %matVLST0
        addVLST0BLADE = (matVLST0(idxVLSTBLADE,:) - matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:)) ...
            * dcmBLADE + matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
        matVLST0 = [matVLST0;addVLST0BLADE];
        
        %matNTVLST0
        addNTVLST0BLADE = (matNTVLST0(idxVLSTBLADE,:) - matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:)) ...
            * dcmBLADE + matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
        matNTVLST0 = [matNTVLST0;addNTVLST0BLADE];
        
        %matCENTER0
        addCENTER0BLADE = (matCENTER0(idxDVEBLADE,:) - matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:)) ...
            * dcmBLADE + matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
        matCENTER0 = [matCENTER0;addCENTER0BLADE];
        
        %valNELE
        valNELE = valNELE + length(idxDVEBLADE);
        
        vecDVEVEHICLE = [vecDVEVEHICLE; vecDVEVEHICLE(idxDVEBLADE)];
        vecDVEWING = [vecDVEWING; vecDVEWING(idxDVEBLADE)];
        vecDVEROTOR = [vecDVEROTOR; vecDVEROTOR(idxDVEBLADE)];
        matSURFACETYPE = [matSURFACETYPE;0,max(matSURFACETYPE(:,2))+1];
        vecDVESURFACE = [vecDVESURFACE;ones(length(idxDVEBLADE),1).*length(matSURFACETYPE(:,1))];
        valPANELS = max(vecDVEPANEL) + 1;
        vecDVEPANEL = [vecDVEPANEL;ones(length(idxDVEBLADE),1).*valPANELS];
        vecDVETIP = [vecDVETIP; vecDVETIP(idxDVEBLADE)];
        vecDVELE = [vecDVELE; vecDVELE(idxDVEBLADE)];
        vecDVETE = [vecDVETE; vecDVETE(idxDVEBLADE)];
        vecDVEROTORBLADE = [vecDVEROTORBLADE; vecDVEROTORBLADE(idxDVEBLADE)];
        vecDVESYM = [vecDVESYM; vecDVESYM(idxDVEBLADE)];
    end
    
    % rotate rotor to axis
    idxVLSTROTOR = unique(matDVE(vecDVEROTOR==n,:));
    idxDVEROTOR = vecDVEROTOR==n; 
    
    matVLST0(idxVLSTROTOR,:)   = matVLST0(idxVLSTROTOR,:)  - matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:);
    matNTVLST0(idxVLSTROTOR,:) = matNTVLST0(idxVLSTROTOR,:)- matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:);
    matCENTER0(idxDVEROTOR,:)  = matCENTER0(idxDVEROTOR,:) - matROTORHUB(n,:) - matVEHORIG(vecROTORVEH(n),:);
    
    % transform rotor from xy plane to hub plane
    dcmROTORHUB = quat2dcm(axang2quat(vrrotvec([0 0 1], matROTORAXIS(n,:))));
    matVLST0(idxVLSTROTOR,:)   = matVLST0(idxVLSTROTOR,:)  * dcmROTORHUB;
    matNTVLST0(idxVLSTROTOR,:) = matNTVLST0(idxVLSTROTOR,:)* dcmROTORHUB;
    matCENTER0(idxDVEROTOR,:)  = matCENTER0(idxDVEROTOR,:) * dcmROTORHUB;
    
    matVLST0(idxVLSTROTOR,:)   = matVLST0(idxVLSTROTOR,:) + matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
    matNTVLST0(idxVLSTROTOR,:) = matNTVLST0(idxVLSTROTOR,:)+ matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
    matCENTER0(idxDVEROTOR,:)  = matCENTER0(idxDVEROTOR,:) + matROTORHUB(n,:) + matVEHORIG(vecROTORVEH(n),:);
    
end






end

