function [matE] = fcnEWING(valNELE, matADJE, vecDVEHVCRD, vecDVELE, vecDVETE)
% OPERA LITE - CHORDWISE CIRCULATION AND VORTICITY
% POSITIVE XSI DIRECTION IS REARWARD, ALIGNED WITH GLOBAL X

%% Finding neighbouring DVEs
idx1 = matADJE(:,2) == 1; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 3; % Borders on 4th local edge (left side)
len = length(idx1(idx1>0,1));

[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);

if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else; idx5 = idx4(eqn_num); % Corresponding edge 3 DVE for every edge 1 DVE
end

%% Vorticity equations between DVEs
cols_zeros = zeros(len,3);

dgamma1 = [cols_zeros ones(len,1) -2.*vecDVEHVCRD(idx3)];
dgamma2 = [cols_zeros ones(len,1) 2.*vecDVEHVCRD(idx5)].*-1;

vort = sparse(len,valNELE*5);
vort = fcnCREATEDSECT(vort, len, 5, idx3, idx5, dgamma1, dgamma2);

%% Circulation equations between DVEs

gamma1 = [ones(len,1) zeros(len,2) -vecDVEHVCRD(idx3) vecDVEHVCRD(idx3).^2];
gamma2 = [ones(len,1) zeros(len,2) vecDVEHVCRD(idx5) vecDVEHVCRD(idx5).^2].*-1;

circ = sparse(len,valNELE*5);
circ = fcnCREATEDSECT(circ, len, 5, idx3, idx5, gamma1, gamma2);

%% Vorticity equations at trailing edges

[edge1,~] = find(vecDVETE == 1);
[edge3,~] = find(vecDVETE == 3);

if isempty(edge1); edge1 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(edge3); edge3 = double.empty(0,1);
end

dgamma1t = [zeros(length(edge1),3) ones(length(edge1),1) -2.*vecDVEHVCRD(edge1)];
dgamma2t = [zeros(length(edge3),3) ones(length(edge3),1) 2.*vecDVEHVCRD(edge3)];
dgammat = [dgamma1t; dgamma2t];

vort_te = sparse(length(edge1) + length(edge3), valNELE*5);
vort_te = fcnCREATEDSECT(vort_te, length(edge1) + length(edge3), 5, [edge1; edge3], [], dgammat, []);

%% Combining all elements into matE
matE = [vort; circ; vort_te];

end

