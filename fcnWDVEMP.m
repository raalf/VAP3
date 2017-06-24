function [ matWDVEMP, matWDVEMPIDX, vecWMPUP, vecWMPDN ] = fcnWDVEMP (matWDVE, matWVLST, matWADJE, valWNELE, vecWDVESYM, vecWDVETIP)
% fcnWPDVEMP takes the wake geometry and works out the shared mid-points
% for relax-wake calculation

% outputs:
% matWDVEMP - ? x 3 list of unique mid-points of wake DVE left/right edge
% matWDVEMPIDX - valWNELE x 2 indeces of mid-points of wake DVE
% vecWMPUP - valWMPNELE x 1 mid-point indeces of upstream mid-point
% vecWMPDN - valWMPNELE x 1 mid-point indeces of downstream mid-point

% Pre-allocation
matWDVEEGL = zeros(valWNELE,2);
matWDVEEGR = zeros(valWNELE,2);

% Find tip, tip right edge equals its own right edge, tip left edge equals its
% own left edge
matWDVEEGR(vecWDVETIP==2,1:2) = [find(vecWDVETIP==2),2.*ones(sum(vecWDVETIP==2),1)];
matWDVEEGL(vecWDVETIP==4,1:2) = [find(vecWDVETIP==4),ones(sum(vecWDVETIP==4),1)];

% find symmetry, symm left edge equals its own left edge, symm right edge equals
% its own right edge
matWDVEEGL(vecWDVESYM==4,1:2) = [find(vecWDVESYM==4),ones(sum(vecWDVESYM==4),1)];
matWDVEEGR(vecWDVESYM==2,1:2) = [find(vecWDVESYM==2),2.*ones(sum(vecWDVESYM==2),1)];

% Find edge with more than 2 panels/
jointADJT = matWADJE(matWADJE(:,4)>0&matWADJE(:,2)==2,:);
matWDVEEGL(jointADJT(:,3),1:2) = [jointADJT(:,1),2.*ones(length(jointADJT(:,1)),1)];

% Rest of the DVEs equals to their own left and right edges
matWDVEEGL(matWDVEEGL(:,2)==0,1:2) = [find(matWDVEEGL(:,1)==0),ones(sum(matWDVEEGL(:,1)==0),1)];
matWDVEEGR(matWDVEEGR(:,2)==0,1:2) = [find(matWDVEEGR(:,1)==0),2*ones(sum(matWDVEEGR(:,1)==0),1)];

% Compare left and right side
matWDVEEG = [matWDVEEGL;matWDVEEGR];

% Lookup DVE coordintes and take the average of 
% points 1,4 for left edges  (matWDVEEG(:,2)==1)
% points 2,3 for right edges (matWDVEEG(:,2)==2)
matWDVEMP = nan(length(matWDVEEG(:,1)),3);
matWDVEMP(matWDVEEG(:,2)==1,1:3) = reshape(mean(reshape(matWVLST(matWDVE(matWDVEEG(matWDVEEG(:,2)==1,1),[1,4])',:)',3,2,[]),2),3,[],1)';
matWDVEMP(matWDVEEG(:,2)==2,1:3) = reshape(mean(reshape(matWVLST(matWDVE(matWDVEEG(matWDVEEG(:,2)==2,1),[2,3])',:)',3,2,[]),2),3,[],1)';

% Filter and index duplicate points to speed up fcnINDVEL
[matWDVEMP,~,C] = unique(matWDVEMP,'rows');
matWDVEMPIDX = reshape(C,valWNELE,2);

valWMPNELE = length(matWDVEMP(:,1));
vecWMPUP = nan(valWMPNELE,1);
vecWMPDN = nan(valWMPNELE,1);


%% Grab upstream Mid-Point Index
% Find upstream DVE element index
dveidxup1 = matWADJE(matWADJE(:,2)==1,1);
dveidxup2 = matWADJE(matWADJE(:,2)==1,3);
mpidxup1 = matWDVEMPIDX(dveidxup1,:);
mpidxup2 = matWDVEMPIDX(dveidxup2,:);
vecWMPUP(mpidxup1(:),1) = mpidxup2(:);
% Find downstream DVE element index
dveidxdn1 = matWADJE(matWADJE(:,2)==3,1);
dveidxdn2 = matWADJE(matWADJE(:,2)==3,3);
mpidxdn1 = matWDVEMPIDX(dveidxdn1,:);
mpidxdn2 = matWDVEMPIDX(dveidxdn2,:);
vecWMPDN(mpidxdn1(:),1) = mpidxdn2(:);





end