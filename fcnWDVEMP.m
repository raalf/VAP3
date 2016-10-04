function [ matWDVEMP, matWDVEMPIDX ] = fcnWDVEMP (matWDVE, matWVLST, matWADJE, valWNELE, vecWDVESYM, vecWDVETIP)
% fcnWPDVEMP takes the wake geometry and works out the shared mid-points
% for relax-wake calculation

% Pre-allocation
matWDVEEGL = zeros(valWNELE,1);
matWDVEEGR = zeros(valWNELE,1);

% Find tip, tip right edge equals its own right edge
% tipRightMP = reshape(mean(reshape(matWVLST(matWDVE(vecWDVETIP==2,[2,3])',:)',3,2,[]),2),3,[],1)';
matWDVEEGR(vecWDVETIP==2,1:2) = [find(vecWDVETIP==2),2.*ones(sum(vecWDVETIP==2),1)]

% find symmetry, symm left edge equals its own left edge
% symLeftMP = reshape(mean(reshape(matWVLST(matWDVE(vecWDVESYM==4,[1,4])',:)',3,2,[]),2),3,[],1)';
matWDVEEGL(vecWDVESYM==4,1:2) = [find(vecWDVESYM==4),ones(sum(vecWDVESYM==4),1)];

% Find edge with more than 2 panels/
jointADJT = matWADJE(matWADJE(:,4)>1&matWADJE(:,2)==2,:);
matWDVEEGL(jointADJT(:,3),1:2) = [jointADJT(:,1),2.*ones(length(jointADJT(:,1)),1)];

% Rest of the DVEs equals to their own left and right edges
matWDVEEGL(matWDVEEGL(:,2)==0,1:2) = [find(matWDVEEGL(:,1)==0),ones(sum(matWDVEEGL(:,1)==0),1)];
matWDVEEGR(matWDVEEGR(:,2)==0,1:2) = [find(matWDVEEGR(:,1)==0),2*ones(sum(matWDVEEGR(:,1)==0),1)];

% Compare left and right side
matWDVEEG = [matWDVEEGL;matWDVEEGR];

% Lookup DVE coordintes and take the average of points 1,4
% for left edges and points 2,3 for right edges
matWDVEMP = nan(length(matWDVEEG(:,1)),3);
matWDVEMP(matWDVEEG(:,2)==1,1:3) = reshape(mean(reshape(matWVLST(matWDVE(matWDVEEG(matWDVEEG(:,2)==1,1),[1,4])',:)',3,2,[]),2),3,[],1)';
matWDVEMP(matWDVEEG(:,2)==2,1:3) = reshape(mean(reshape(matWVLST(matWDVE(matWDVEEG(matWDVEEG(:,2)==2,1),[2,3])',:)',3,2,[]),2),3,[],1)';

% Filter and index duplicate points to speed up fcnINDVEL
[matWDVEMP,~,C] = unique(matWDVEMP,'rows');
matWDVEMPIDX = reshape(C,valWNELE,2);

end