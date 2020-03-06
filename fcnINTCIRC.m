function [intcirc] = fcnINTCIRC(SURF)

idx1 = SURF.vecDVELE == 1; %index of LE vectors (will be the same)

A = zeros(1,SURF.valNELE);
B = zeros(1,SURF.valNELE);
C = zeros(1,SURF.valNELE);

A(idx1) = SURF.matCOEFF(idx1,1);
B(idx1) = SURF.matCOEFF(idx1,2);
C(idx1) = SURF.matCOEFF(idx1,3);
% if any other row, A= A-Aupstream, B= B-Bupstream, C= C-Cupstream

idx2 = SURF.vecDVELE == 1; %idx2 since we need to do this even for triangles
dvenum = find(idx2==0); %dvenum in question
idxf = SURF.matADJE((ismember(SURF.matADJE(:,1), dvenum) & SURF.matADJE(:,2) == 1),3); %upstream dve num
% A(idx2 ==0) = (SURF.matCOEFF(idx2==0,1)-SURF.matCOEFF(idxf,1));
% B(idx2 ==0) = (SURF.matCOEFF(idx2==0,2)-SURF.matCOEFF(idxf,2));
% C(idx2 ==0) = (SURF.matCOEFF(idx2==0,3)-SURF.matCOEFF(idxf,3));

% intcirc = sum(((A .*2 .* SURF.vecDVEHVSPN'+  C./3.*2.*SURF.vecDVEHVSPN'.*SURF.vecDVEHVSPN'.*SURF.vecDVEHVSPN'))',1); % Integrated circulation across DVE

A = SURF.matCOEFF(:,1);
C = SURF.matCOEFF(:,3);

intcirc = sum(((A .*2 .* SURF.vecDVEHVSPN + C./3.*2.*SURF.vecDVEHVSPN.*SURF.vecDVEHVSPN.*SURF.vecDVEHVSPN)),1); % Integrated circulation across DVE

end