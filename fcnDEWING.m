function matD = fcnDEWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP, vecN, vecDVEHVCRD, vecDVELE, vecDVETE)

%% Finding neighbouring DVEs (Spanwise)
idx1 = matADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 4; % Borders on 4th local edge (left side)
len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);

if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else; idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE
end

%% Vorticity equations between DVEs in Spanwise Direction
% dgamma1 = [2.*vecDVEHVSPN(idx3) ones(len,1) zeros(len,4); ... % Spanwise
%           zeros(len,4), ones(len,1) zeros(len,1)];                          % Chordwise, at midchord    
% dgamma2 = [-2.*vecDVEHVSPN(idx5) ones(len,1) zeros(len,4); ... % Spanwise
%           zeros(len,4), ones(len,1) zeros(len,1)].*-1;
% vort_span = fcnCREATEDSECT(sparse(len*2,valNELE*6), len*2, 6, [idx3; idx3], [idx5; idx5], dgamma1, dgamma2);

% dgamma1 = [2.*vecDVEHVSPN(idx3) ones(len,1) zeros(len,2) ones(len,1) zeros(len,1)];   
% dgamma2 = [-2.*vecDVEHVSPN(idx5) ones(len,1) zeros(len,2) ones(len,1) zeros(len,1)].*-1;
%     
% vort_span = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, dgamma1, dgamma2);

dgamma1 = [2.*vecDVEHVSPN(idx3) ones(len,1) zeros(len,4)];    
dgamma2 = [-2.*vecDVEHVSPN(idx5) ones(len,1) zeros(len,4)].*-1;

vort_span = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, dgamma1, dgamma2);

%% Circulation equations between DVEs in Spanwise Direction
gamma1 = [vecDVEHVSPN(idx3).^2 vecDVEHVSPN(idx3) ones(len,1) zeros(len,2) ones(len,1)];
gamma2 = [vecDVEHVSPN(idx5).^2 -vecDVEHVSPN(idx5) ones(len,1) zeros(len,2) ones(len,1)].*-1;

gamma_span = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, gamma1, gamma2);      

%% Finding neighbouring DVEs (Chordwise)
idx1 = matADJE(:,2) == 1; % Borders on 1st local edge (right side)
idx2 = matADJE(:,2) == 3; % Borders on 3rd local edge (left side)
len = length(idx1(idx1>0,1));

[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);

if isempty(idx3); idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else; idx5 = idx4(eqn_num); % Corresponding edge 3 DVE for every edge 1 DVE
end

%% Vorticity equations between DVEs in Chordwise Direction
% dgamma1 = [zeros(len,1) ones(len,1) zeros(len,4); ... % Spanwise, at midspan
%           zeros(len,3) -2.*vecDVEHVCRD(idx3) ones(len,1) zeros(len,1)]; % Chordwise   
% dgamma2 = [zeros(len,1) ones(len,1) zeros(len,4); ... % Spanwise, at midspan
%           zeros(len,3) 2.*vecDVEHVCRD(idx5) ones(len,1) zeros(len,1)].*-1; % Chordwise      
% vort_chord = fcnCREATEDSECT(sparse(len*2,valNELE*6), len*2, 6, [idx3; idx3], [idx5; idx5], dgamma1, dgamma2);

% dgamma1 = [zeros(len,1) ones(len,1) zeros(len,1) -2.*vecDVEHVCRD(idx3) ones(len,1) zeros(len,1)];
% dgamma2 = [zeros(len,1) ones(len,1) zeros(len,1) 2.*vecDVEHVCRD(idx5) ones(len,1) zeros(len,1)].*-1;
%       
% vort_chord = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, dgamma1, dgamma2);

dgamma1 = [zeros(len,3) -2.*vecDVEHVCRD(idx3) ones(len,1) zeros(len,1)]; % Chordwise      
dgamma2 = [zeros(len,3) 2.*vecDVEHVCRD(idx5) ones(len,1) zeros(len,1)].*-1; % Chordwise

vort_chord = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, dgamma1, dgamma2);

%% Circulation equations between DVEs in Chordwise Direction
gamma1 = [zeros(len,2) ones(len,1) vecDVEHVCRD(idx3).^2 -vecDVEHVCRD(idx3) ones(len,1)];
gamma2 = [zeros(len,2) ones(len,1) vecDVEHVCRD(idx5).^2 vecDVEHVCRD(idx5) ones(len,1)].*-1;

gamma_chord = fcnCREATEDSECT(sparse(len,valNELE*6), len, 6, idx3, idx5, gamma1, gamma2); 

%% Circulation equations at wingtip
[tip2,~] = find(vecDVETIP == 2);
[tip4,~] = find(vecDVETIP == 4);

if vecN == 1; tip2 = [1:length(vecDVETIP)]'; tip4 = [1:length(vecDVETIP)]'; end

len = length(tip2) + length(tip4);

if isempty(tip2); tip2 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(tip4); tip4 = double.empty(0,1);
end

gamma1t = [vecDVEHVSPN(tip2).^2 vecDVEHVSPN(tip2) ones(length(tip2),1) zeros(length(tip2),2) ones(length(tip2),1)];
gamma2t = [vecDVEHVSPN(tip4).^2 -vecDVEHVSPN(tip4) ones(length(tip4),1) zeros(length(tip4),2) ones(length(tip4),1)];
gammat = [gamma1t; gamma2t];

circ_tip = fcnCREATEDSECT(sparse(len, valNELE*6), len, 6, [tip2; tip4], [], gammat, []);

%% Vorticity equations at trailing edges
[edge1,~] = find(vecDVETE == 1);
[edge3,~] = find(vecDVETE == 3);

if isempty(edge1); edge1 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(edge3); edge3 = double.empty(0,1);
end

dgamma1t = [zeros(length(edge1),3) -2.*vecDVEHVCRD(edge1) ones(length(edge1),1) zeros(length(edge1),1)];
dgamma2t = [zeros(length(edge3),3) 2.*vecDVEHVCRD(edge3) ones(length(edge3),1) zeros(length(edge3),1)];
dgammat = [dgamma1t; dgamma2t];

vort_te = fcnCREATEDSECT(sparse(length(edge1) + length(edge3), valNELE*6), length(edge1) + length(edge3), 6, [edge1; edge3], [], dgammat, []);

%% Combining all together
matD = [gamma_span; gamma_chord; vort_span; vort_chord; circ_tip; vort_te];

end