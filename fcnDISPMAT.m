function [ matVUP, matVNOW, matVDOWN] = fcnDISPMAT(matWDVEMPIND, vecWMPUP, vecWMPDN );
%%
valWMPNELE = length(matWDVEMPIND(:,1));
matVUP = nan(valWMPNELE,3);
matVNOW = nan(valWMPNELE,3);
matVDOWN = nan(valWMPNELE,3);
% case 1 - Lastest Wake Elements, no upstream elements
idx1 = isnan(vecWMPUP);
% (Current, [0,0,0], Downstream)
matVUP(idx1,1:3) = matWDVEMPIND(idx1,1:3);
matVNOW(idx1,1:3) = zeros(sum(idx1),3);
matVDOWN(idx1,1:3) = matWDVEMPIND(vecWMPDN(idx1),1:3);

% case 2 - Oldest Wake Elements, DO NOT RELAX
% ([0,0,0],[0,0,0],[0,0,0])
idx2 = isnan(vecWMPDN);
% matVUP(idx2,1:3) = zeros(sum(idx2),3);
% matVNOW(idx2,1:3) = zeros(sum(idx2),3);
% matVDOWN(idx2,1:3) = zeros(sum(idx2),3);
matVUP(idx2,1:3) = matWDVEMPIND(vecWMPUP(idx2),1:3);
matVNOW(idx2,1:3) = matWDVEMPIND(vecWMPUP(idx2),1:3);
matVDOWN(idx2,1:3) = matWDVEMPIND(vecWMPUP(idx2),1:3);

% case 3 - 2nd Oldest Wake Elements
% idx3, find downstream MP which is on the idx2 list
% (Current, Current, Current)
idx3 = ismember(vecWMPDN,find(idx2));
matVUP(idx3,1:3) = matWDVEMPIND(idx3,1:3);
matVNOW(idx3,1:3) = matWDVEMPIND(idx3,1:3);
matVDOWN(idx3,1:3) = matWDVEMPIND(idx3,1:3);

% case 4 - the rest of the MP list
% idx4, MP not in any idx list before
idx4 = ~idx1&~idx2&~idx3;
% (Upstream, Current, Downstream)
matVUP(idx4,1:3) = matWDVEMPIND(vecWMPUP(idx4),1:3);
matVNOW(idx4,1:3) = matWDVEMPIND(idx4,1:3);
matVDOWN(idx4,1:3) = matWDVEMPIND(vecWMPDN(idx4),1:3);

if sum(sum(isnan([matVUP;matVNOW;matVDOWN]))) > 0
    warning('fcnDISPMAT: Displace Matrix Error');
end

end