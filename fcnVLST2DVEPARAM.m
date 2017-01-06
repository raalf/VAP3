function [ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, matVLST, matCENTER ] = fcnVLST2DVEPARAM( matDVE, matNPVLST )

%FCNVLST2DVEPARAMS Summary of this function goes here
%   Detailed explanation goes here

P1 = matNPVLST(matDVE(:,1),:);
P2 = matNPVLST(matDVE(:,2),:);
P3 = matNPVLST(matDVE(:,3),:);
P4 = matNPVLST(matDVE(:,4),:);

matCENTER = (P1+P2+P3+P4)/4;


[ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, ...
    matVLST, ~, ~, ~] = fcnDVECORNER2PARAM( matCENTER, P1, P2, P3, P4 );



end

