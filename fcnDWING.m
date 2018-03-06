function [matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP, vecN)
% Currently creates the upper 2/3rds of the D-matrix, which are the boundary conditions based
% on adjacent DVEs, etc.
%
% INPUT:
%   valNELE - total number of DVEs
%   matADJE - ? x 3 adjacency matrix, where columns are: DVE | local edge | adjacent DVE
%   vecDVEHVSPN - valNELE x 1 vector of DVE half spans
%   vecDVESYM - valNELE x 1 vector of which DVEs have symmetry on which edge (0 for no symmetry,
%                   2 for local edge 2, 4 for local edge 4)
%   vecDVETIP - valNELE x 1 vector of which DVEs are at the wingtip. Similar format to vecDVESYM above
% OUTPUT:
%   matD - currently, upper 2/3 of D-matrix

%% Finding neighbouring DVEs
idx1 = matADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 4; % Borders on 4th local edge (left side)
len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
%[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
%eqn_num = cell2mat(LB);
[~, eqn_num] = ismember([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'rows'); % this line is added to replace ismembertol to suppport UINT8 variable type

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);

if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else; idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE
end

%% Vorticity equations between DVEs
dgamma1 = [zeros(len,1) ones(len,1) 2.*vecDVEHVSPN(idx3)];
dgamma2 = [zeros(len,1) ones(len,1) -2.*vecDVEHVSPN(idx5)].*-1;

vort = sparse(len,valNELE*3);
vort = fcnCREATEDSECT(vort, len, 3, idx3, idx5, dgamma1, dgamma2);

%% Circulation equations between DVEs
gamma1 = [ones(len,1) vecDVEHVSPN(idx3) vecDVEHVSPN(idx3).^2];
gamma2 = [ones(len,1) -vecDVEHVSPN(idx5) vecDVEHVSPN(idx5).^2].*-1;

circ = sparse(len,valNELE*3);
circ = fcnCREATEDSECT(circ, len, 3, idx3, idx5, gamma1, gamma2);

%% Vorticity equations at symmetry
vort_sym = [];
if any(vecDVESYM)
    [idx10, ~] = find(vecDVESYM);
    len = length(idx10);
    locedge = vecDVESYM(idx10);
    
    dgamma_sym = zeros(len,3);
    dgamma_sym(locedge == 2,:) = [zeros(length(idx10(locedge == 2)),1) ones(length(idx10(locedge == 2)),1) 2.*vecDVEHVSPN(idx10(locedge == 2))];
    dgamma_sym(locedge == 4,:) = [zeros(length(idx10(locedge == 4)),1) ones(length(idx10(locedge == 4)),1) -2.*vecDVEHVSPN(idx10(locedge == 4))];
    
    vort_sym = sparse(len,valNELE*3);
    vort_sym = fcnCREATEDSECT(vort_sym, len, 3, uint32(idx10), [], dgamma_sym, []);
end

%% Circulation equations at wingtip
[tip2,~] = find(vecDVETIP == 2);
[tip4,~] = find(vecDVETIP == 4);

tip2 = uint32(tip2);
tip4 = uint32(tip4);

if sum(vecN) == 1 && ~any(vecDVESYM); tip2 = [1:length(vecDVETIP)]'; tip4 = [1:length(vecDVETIP)]'; end

len = length(tip2) + length(tip4);

if isempty(tip2); tip2 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(tip4); tip4 = double.empty(0,1);
end

gamma1t = [ones(length(tip2),1) vecDVEHVSPN(tip2) vecDVEHVSPN(tip2).^2];
gamma2t = [ones(length(tip4),1) -vecDVEHVSPN(tip4) vecDVEHVSPN(tip4).^2];
gammat = [gamma1t; gamma2t];

circ_tip = sparse(len, valNELE*3);
circ_tip = fcnCREATEDSECT(circ_tip, len, 3, [tip2; tip4], [], gammat, []);

%% Combining all elements into matD
matD = [vort; circ; vort_sym; circ_tip];

end