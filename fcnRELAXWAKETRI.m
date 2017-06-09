function [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
    matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKETRI(vecUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
    matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
    vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING, vecDVETE, flagTRI)


%% Moving vertices

matWDVEMP = matWVLST;

% Get induced velocity at vertices
[ matWDVEMPIND ] = fcnINDVEL(matWDVEMP,valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, flagTRI);


% indicies to not move (attached to wing)
 newest_row = sort([valWNELE:-2:valWNELE-valWSIZE*2+1]-1)';
 idx1 = unique(reshape(matWDVE(newest_row,1:2),[],1),'rows');
 
 oldest_row = [2:2:valWSIZE*2];
 idx2 = unique(reshape(matWDVE(oldest_row,1:4),[],1),'rows');
 
 idx3 = ones(size(matWDVEMP,1),1);
idx3(idx1)= 0;
idx3(idx2) = 0;
%Move verticies by indvel*time
matWVLST(idx3~=0,:) = matWDVEMP(idx3~=0,:) + matWDVEMPIND(idx3~=0,:).*valDELTIME;

%find new collocaiton points
WP1 = matWVLST(matWDVE(:,1),:);
WP2 = matWVLST(matWDVE(:,2),:);
WP3 = matWVLST(matWDVE(:,3),:);
WP4 = matWVLST(matWDVE(:,4),:);
matWCENTER = (WP1(:,:)+WP2(:,:)+WP3(:,:)+WP4(:,:))./4;

% update relax wake dves
[vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWDVENORM, ...
    matWVLST, matWDVE, ~, idxWVLST] = fcnDVECORNER2PARAM( matWCENTER, WP1, WP2, WP3, WP4 );

% For singularity factor updating, each row of wake elements needs to be seen as its own "wing"
% so we are adding on to the vecWDVEWING by timestep number to create the right form to pass into
% the surface singularity factor function
tswing = vecWDVEWING + reshape(repmat([0:max(vecWDVEWING):((valWNELE/(valWSIZE*2))-1)*max(vecWDVEWING)], valWSIZE*2,1),[],1);

% Updating wake singularity factor
[vecWK] = fcnSINGFCT(valWNELE, tswing, vecWDVETIP, vecWDVEHVSPN);

end

