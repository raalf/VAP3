function [matD] = fcnDWING(SURF, INPU)
% Currently creates the upper 2/3rds of the D-matrix, which are the boundary conditions based
% on adjacent DVEs, etc.
%
% INPUT:
%   SURF.valNELE - total number of DVEs
%   SURF.matADJE - ? x 3 adjacency matrix, where columns are: DVE | local edge | adjacent DVE
%   SURF.vecDVEHVSPN - SURF.valNELE x 1 vector of DVE half spans
%   SURF.vecDVESYM - SURF.valNELE x 1 vector of which DVEs have symmetry on which edge (0 for no symmetry,
%                   2 for local edge 2, 4 for local edge 4)
%   SURF.vecDVETIP - SURF.valNELE x 1 vector of which DVEs are at the wingtip. Similar format to SURF.vecDVESYM above
% OUTPUT:
%   matD - currently, upper 2/3 of D-matrix

%% Finding neighbouring DVEs (not splits)
idx1 = SURF.matADJE(:,2) == 2 & SURF.matADJE(:,4) == 1; % Borders on 2nd local edge (right side)
idx2 = SURF.matADJE(:,2) == 4 & SURF.matADJE(:,4) == 1; % Borders on 4th local edge (left side)
len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
%[~, LB] = ismembertol([SURF.matADJE(idx1,1) SURF.matADJE(idx1,3)],[SURF.matADJE(idx2,3) SURF.matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
%eqn_num = cell2mat(LB);
[~, eqn_num] = ismember([SURF.matADJE(idx1,1) SURF.matADJE(idx1,3)],[SURF.matADJE(idx2,3) SURF.matADJE(idx2,1)],'rows'); % this line is added to replace ismembertol to suppport UINT8 variable type

idx3 = SURF.matADJE(idx1,1);
idx4 = SURF.matADJE(idx2,1);

if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else; idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE
end

%% Vorticity equations between DVEs
dgamma1 = [zeros(len,1) ones(len,1) 2.*SURF.vecDVEHVSPN(idx3)];
dgamma2 = [zeros(len,1) ones(len,1) -2.*SURF.vecDVEHVSPN(idx5)].*-1;

vort = sparse(len,SURF.valNELE*3);
vort = fcnCREATEDSECT(vort, len, 3, idx3, idx5, [], dgamma1, dgamma2, []);

%% Circulation equations between DVEs
gamma1 = [ones(len,1) SURF.vecDVEHVSPN(idx3) SURF.vecDVEHVSPN(idx3).^2];
gamma2 = [ones(len,1) -SURF.vecDVEHVSPN(idx5) SURF.vecDVEHVSPN(idx5).^2].*-1;

circ = sparse(len,SURF.valNELE*3);
circ = fcnCREATEDSECT(circ, len, 3, idx3, idx5, [], gamma1, gamma2, []);

%% Splits
vort_split = [];
circ_split = [];
idx1 = SURF.matADJE(:,2) == 2 & SURF.matADJE(:,4) == 2; % Borders on 2nd local edge (right side)
idx2 = SURF.matADJE(:,2) == 4 & SURF.matADJE(:,4) == 2; % Borders on 4th local edge (left side)

if any(idx1)
    % Mapping equations evaluated on edge 2 to equations evaluated on edge 4
    [~, eqn_num] = ismember([SURF.matADJE(idx1,1) SURF.matADJE(idx1,3)],[SURF.matADJE(idx2,3) SURF.matADJE(idx2,1)],'rows'); % this line is added to replace ismembertol to suppport UINT8 variable type

    idx3 = SURF.matADJE(idx1,1);
    idx4 = SURF.matADJE(idx2,1);
    if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
    else; idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE
    end

    % Vorticity
    len = length(idx5);
    dgamma1s = [zeros(len,1) ones(len,1) 2.*SURF.vecDVEHVSPN(idx3)];
    dgamma2s = [zeros(len,1) ones(len,1) -2.*SURF.vecDVEHVSPN(idx5)].*-1;
    vort_split = fcnCREATEDSECT(sparse(len,SURF.valNELE*3), len, 3, idx3, idx5, [], dgamma1s, dgamma2s, []);
    
    % Circulation
    ldve = reshape(unique(idx3),1,1,[]);

    tmp = repmat(idx3,1,1,length(ldve)) == repmat(ldve,length(idx5),1,1);
    tmp2 = repmat(idx5,1,1,length(ldve));
    rdve = reshape(tmp2(tmp),2,1,[]);
    
    gamma_ldve = [ones(size(ldve,3),1) SURF.vecDVEHVSPN(ldve) SURF.vecDVEHVSPN(ldve).^2];
    gamma_rdve = [ones(size(rdve,1),1,size(ldve,3)) SURF.vecDVEHVSPN(rdve) SURF.vecDVEHVSPN(rdve).^2].*-1;
    gamma_rdve = permute(gamma_rdve, [3 2 1]); 

    len = size(ldve,3); 
    ldve = reshape(ldve,[],1,1);
    rdve = permute(rdve, [3 2 1]);

    circ_split = fcnCREATEDSECT(sparse(len, SURF.valNELE*3), len, 3, ldve, rdve(:,:,1), rdve(:,:,2), gamma_ldve, gamma_rdve(:,:,1), gamma_rdve(:,:,2)); 
end

%% Vorticity equations at symmetry
vort_sym = [];
if any(SURF.vecDVESYMEDGE)
    [idx10, ~] = find(SURF.vecDVESYMEDGE);
    len = length(idx10);
    locedge = SURF.vecDVESYMEDGE(idx10);
    
    dgamma_sym = zeros(len,3);
    dgamma_sym(locedge == 2,:) = [zeros(length(idx10(locedge == 2)),1) ones(length(idx10(locedge == 2)),1) 2.*SURF.vecDVEHVSPN(idx10(locedge == 2))];
    dgamma_sym(locedge == 4,:) = [zeros(length(idx10(locedge == 4)),1) ones(length(idx10(locedge == 4)),1) -2.*SURF.vecDVEHVSPN(idx10(locedge == 4))];
    
    vort_sym = sparse(len,SURF.valNELE*3);
    vort_sym = fcnCREATEDSECT(vort_sym, len, 3, uint32(idx10), [], [], dgamma_sym, [], []);
end

%% Circulation equations at wingtip
[tip2,~] = find(SURF.vecDVETIP == 2);
[tip4,~] = find(SURF.vecDVETIP == 4);

tip2 = uint32(tip2);
tip4 = uint32(tip4);

if sum(INPU.vecN) == 1 && ~any(SURF.vecDVESYMEDGE); tip2 = [1:length(SURF.vecDVETIP)]'; tip4 = [1:length(SURF.vecDVETIP)]'; end

len = length(tip2) + length(tip4);

if isempty(tip2); tip2 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(tip4); tip4 = double.empty(0,1);
end

gamma1t = [ones(length(tip2),1) SURF.vecDVEHVSPN(tip2) SURF.vecDVEHVSPN(tip2).^2];
gamma2t = [ones(length(tip4),1) -SURF.vecDVEHVSPN(tip4) SURF.vecDVEHVSPN(tip4).^2];
gammat = [gamma1t; gamma2t];

circ_tip = sparse(len, SURF.valNELE*3);
circ_tip = fcnCREATEDSECT(circ_tip, len, 3, [tip2; tip4], [], [], gammat, [], []);

%% Combining all elements into matD
matD = [vort; circ; vort_sym; vort_split; circ_split; circ_tip];

end