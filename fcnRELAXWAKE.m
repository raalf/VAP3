function [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
    matWVLST, matWDVE, idxWVLST, vecWK] = fcnRELAXWAKE(matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
    matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVEHVSPN, vecDVELESWP, ...
    vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVELESWP, vecWDVEPITCH, ...
    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING)
%FCNRLXWAKE Summary of this function goes here
%   Detailed explanation goes here
[ matWDVEMP, matWDVEMPIDX, vecWMPUP, vecWMPDN ] = fcnWDVEMP(matWDVE, matWVLST, matWADJE, valWNELE, vecWDVESYM, vecWDVETIP);

% Get mid-points induced velocity
[ matWDVEMPIND ] = fcnINDVEL(matWDVEMP,valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);

% Assemble matrices for fcnDISPLACE (vup, vnow, vdown)
[ matVUP, matVNOW, matVDOWN ] = fcnDISPMAT(matWDVEMPIND, vecWMPUP, vecWMPDN );

% Calculate mid-point displacement
[matWDVEMPRLX] = fcnDISPLACE(matVUP, matVNOW, matVDOWN, matWDVEMP, valDELTIME);

% update matWCENTER
matWCENTER = matWDVEMPRLX(matWDVEMPIDX(:,1),:)+(0.5*(matWDVEMPRLX(matWDVEMPIDX(:,2),:)-matWDVEMPRLX(matWDVEMPIDX(:,1),:)));

% Calculate the leading edge mid-points by refering to upstream dves
matWDVELEMPIDX = flipud(reshape(1:valWNELE,valWSIZE,[])');
matWDVELEMP = nan(valWNELE,3);
for wakerow = 1:valTIMESTEP
    if wakerow == 1% || wakerow == valTIMESTEP % freshest row of wake, closest to wing TE
        matWDVELEMP(matWDVELEMPIDX(wakerow,:),(1:3)) = permute(mean(reshape(matWVLST(matWDVE(matWDVELEMPIDX(wakerow,:),[1,2]),:)',3,[],2),3),[3 2 1]);
    else % rest of the wake dves
        prevLE = matWDVELEMP(matWDVELEMPIDX(wakerow-1,:),(1:3));
        prevXO = matWCENTER(matWDVELEMPIDX(wakerow-1,:),(1:3));
        matWDVELEMP(matWDVELEMPIDX(wakerow,:),(1:3)) = prevLE + 2.*(prevXO - prevLE);
    end
end

% half chord vector, leading edge midpoint -> control point
crdvec = matWCENTER - matWDVELEMP;
% half span vector, control point -> right edge midpoint
spnvec = matWDVEMPRLX(matWDVEMPIDX(:,2),:) - matWCENTER;

% fix oldest wake
semiinfvec = matWCENTER(matWDVELEMPIDX(end-1,:),:)-matWDVELEMP(matWDVELEMPIDX(end-1,:),:);
matWCENTER(matWDVELEMPIDX(end,:),:) = matWDVELEMP(matWDVELEMPIDX(end,:),:)+semiinfvec;
crdvec(matWDVELEMPIDX(end,:),:) = crdvec(matWDVELEMPIDX(end-1,:),:);
spnvec(matWDVELEMPIDX(end,:),:) = spnvec(matWDVELEMPIDX(end-1,:),:);

% Calculate four corner point by adding vectors to control point coordinates
WP1 = matWCENTER - crdvec - spnvec;
WP2 = matWCENTER - crdvec + spnvec;
WP3 = matWCENTER + crdvec + spnvec;
WP4 = matWCENTER + crdvec - spnvec;
%%

%%
% update relax wake dves
[vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWDVENORM, ...
    matWVLST, matWDVE, ~, idxWVLST] = fcnDVECORNER2PARAM( matWCENTER, WP1, WP2, WP3, WP4 );

% For singularity factor updating, each row of wake elements needs to be seen as its own "wing"
% so we are adding on to the vecWDVEWING by timestep number to create the right form to pass into
% the surface singularity factor function
tswing = vecWDVEWING + reshape(repmat([0:max(vecWDVEWING):((valWNELE/valWSIZE)-1)*max(vecWDVEWING)], valWSIZE,1),[],1);

% Updating wake singularity factor
[vecWK] = fcnSINGFCT(valWNELE, tswing, vecWDVETIP, vecWDVEHVSPN);

end

