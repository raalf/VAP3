function [ matADJE, vecDVESYM, vecDVETIP, vecDVELE, vecDVETE ] = fcnDVEADJT( imP1, imP2, imP3, imP4, valNELE, vecDVEPANEL, vecSYM )
%FCNDVEADJT Summary of this function goes here
%   Detailed explanation goes here

%Pre-allocation
vecDVESYM   = zeros(valNELE,1);
vecDVETIP   = zeros(valNELE,1);
vecDVELE    = zeros(valNELE,1);
vecDVETE    = zeros(valNELE,1);

[~,~,idxNPV] = unique([imP1;imP2;imP3;imP4],'rows');
idxNPV = reshape(idxNPV,valNELE,4);
temp = sort([idxNPV(:,[1,2]);idxNPV(:,[2,3]);idxNPV(:,[3,4]);idxNPV(:,[4,1])],2); %sort to ensure the edge is align to same direction
[~,~,j2] = unique(temp,'rows','stable');
EIDX = reshape(j2,valNELE,4);

%row: DVE# | Local Edge # | Glob. edge#
j = [repmat((1:valNELE)',4,1),reshape(repmat(1:4,valNELE,1),valNELE*4,1),reshape(EIDX,valNELE*4,1)];
[j1,~] = histc(j(:,3),unique(j(:,3)));
j = [j,j1(j(:,3))-1];
%Currently the procedure was done in two for loops. May be modified in
%later days if performance improvement is required.
matADJE = uint32(nan(sum(j(:,4)),4));
k = j(j(:,4)~=0,:);
c = 0;

for i = 1:length(k(:,1))
    for i2 = 1:k(i,4)
        c = c+1;
        currentdve = k(i,1);
        currentlocaledge = k(i,2);
        currentedge = k(i,3);
        
        dvefulllist = j(j(:,3)==currentedge,1);
        dvefilterlist = dvefulllist(dvefulllist~=currentdve);
        matADJE(c,:) = [currentdve currentlocaledge dvefilterlist(i2) k(i,4)];
    end
end

%sort matADJE by dve#
[~,tempB] = sort(matADJE(:,1));
matADJE = matADJE(tempB,:);


%% Handle symmetry and wing tip DVEs
idx1 = (j(:,4)==0&(j(:,2)==2|j(:,2)==4));
dveedge2panel = repmat(vecDVEPANEL,4,1);

findTIPSYM = [j(idx1,:),dveedge2panel(idx1,:)];
tempTIP = nan(length(findTIPSYM(:,1)),1);
% If the panel has vecSYM = 1, those panels' local edge 4 has symmetry
panelidx = find(vecSYM==1);
symidx = (findTIPSYM(:,2)==4 & ismember(findTIPSYM(:,5),panelidx));
tempTIP(symidx) = 1;
dveidx = findTIPSYM(symidx,1);
vecDVESYM(dveidx) = 4;

% If the panel has vecSYM = 2, those panels' local edge 2 has symmetry
panelidx = find(vecSYM==2);
symidx = (findTIPSYM(:,2)==2 & ismember(findTIPSYM(:,5),panelidx));
tempTIP(symidx) = 1;
dveidx = findTIPSYM(symidx,1);
vecDVESYM(dveidx) = 2;

% If the edge is not touching another dve nor symmetry edge,
% Define it as wing tip
% WHAT IF BOTH EDGES ARE TIPS???????????????????????????????????????????????????
tipidx = isnan(tempTIP);
dveidx=  findTIPSYM(tipidx,1);
vecDVETIP(dveidx) = findTIPSYM(tipidx,2);

% Get DVE index where TE appears if col2=1 & col4=0;
% use matrix j, col1:dve, col2:Local.edge, col3:Glob.edge, col4:
leIdx = j(j(:,2)==1&j(:,4)==0,1);
vecDVELE(leIdx) = 1;

% Get DVE index where TE appears if col2=3 & col4=0;
% use matrix j, col1:dve, col2:Local.edge, col3:Glob.edge, col4:
teIdx = j(j(:,2)==3&j(:,4)==0,1);
vecDVETE(teIdx) = 3;

end

