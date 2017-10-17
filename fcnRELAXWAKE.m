function [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
    matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKE(matUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
    matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVETE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
    vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING, flagSTEADY)
%FCNRLXWAKE Summary of this function goes here
%   Detailed explanation goes here
[ matWDVEMP, matWDVEMPIDX, vecWMPUP, vecWMPDN ] = fcnWDVEMP(matWDVE, matWVLST, matWADJE, valWNELE, vecWDVESYM, vecWDVETIP);

% Get mid-points induced velocity
[ matWDVEMPIND ] = fcnINDVEL(matWDVEMP,valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, 0, flagSTEADY);

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
    if wakerow == 1 %|| wakerow == valTIMESTEP % freshest row of wake, closest to wing TE
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
% matWCENTER(matWDVELEMPIDX(end,:),:) = matWDVELEMP(matWDVELEMPIDX(end,:),:)+semiinfvec;

xsi0 = repmat(sqrt(semiinfvec(:,1).^2+semiinfvec(:,2).^2+semiinfvec(:,3).^2),1,3);
matWCENTER(matWDVELEMPIDX(end,:),:) = matWCENTER(matWDVELEMPIDX(end-1,:),:)+semiinfvec+matUINF(vecDVETE == 3).*xsi0;
% crdvec(matWDVELEMPIDX(end,:),:) = crdvec(matWDVELEMPIDX(end-1,:),:);
% spnvec(matWDVELEMPIDX(end,:),:) = spnvec(matWDVELEMPIDX(end-1,:),:);


% Calculate four corner point by adding vectors to control point coordinates
WP1 = matWCENTER - crdvec - spnvec;
WP2 = matWCENTER - crdvec + spnvec;
WP3 = matWCENTER + crdvec + spnvec;
WP4 = matWCENTER + crdvec - spnvec;

%% Recalculating oldest row 4 corner points
% This overwrites the WP1-WP4 points of oldest wake elements
oldestwake = reshape(matWDVELEMPIDX(end,:),[],1);
secondoldestwake = reshape(matWDVELEMPIDX(end-1,:),[],1);
WP1(oldestwake,:) = WP4(secondoldestwake,:);
WP2(oldestwake,:) = WP3(secondoldestwake,:);
timesteptranslate = matUINF(vecDVETE == 3,:)*valDELTIME;
WP4(oldestwake,:) = WP1(oldestwake,:)+timesteptranslate;
WP3(oldestwake,:) = WP2(oldestwake,:)+timesteptranslate;
matWCENTER(oldestwake,:) = (WP1(oldestwake,:)+WP2(oldestwake,:)+WP3(oldestwake,:)+WP4(oldestwake,:))./4;



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

