function [matWD, vecWR] = fcnWDWAKE(all_DVEs, matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM, vecN)
% Creates the Wake D-matrix, for updating the vorticity coefficients
% for the wake to account for stretching.

% INPUT:
%   dves - list of wake DVE numbers we are creating the matrix for
%   matWADJE - ? x 3 adjacency matrix, where columns are: Wake DVE | local edge | adjacent Wake DVE
%   vecWDVEHVSPN - nelements x 1 vector of wake DVE half spans
%   vecWDVESYM - nelements x 1 vector of which wake DVEs have symmetry on which edge (0 for no symmetry,
%                   2 for local edge 2, 4 for local edge 4)
%   vecWDVETIP - nelements x 1 vector of which wake DVEs are at the wingtip. Similar format to vecWDVESYM above
% OUTPUT:
%   matWD - wake D-matrix

nelements = length(all_DVEs);

%% Vorticity and circulation equations between DVEs ----------------------------------------------------------------------

idx1 = matWADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx2 = matWADJE(:,2) == 4; % Borders on 4th local edge (left side)

len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
% This returns a vector where we can see for edge 2 equations, which is the
% matching edge 4 equation
% This will work well for split vorticity but it will need modification for split circulation later on
% [~, LB] = ismembertol([matWADJE(idx2,3) matWADJE(idx2,1)],[matWADJE(idx1,1) matWADJE(idx1,3)],'ByRows',true,'OutputAllIndices',true);
% [~, LB] = ismembertol([matWADJE(idx1,1) matWADJE(idx1,3)],[matWADJE(idx2,3) matWADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
% eqn_num = cell2mat(LB);
[~, eqn_num] = ismember([matWADJE(idx1,1) matWADJE(idx1,3)],[matWADJE(idx2,3) matWADJE(idx2,1)],'rows'); % this line is added to replace ismembertol to suppport UINT8 variable type

idx3 = matWADJE(idx1,1);
idx4 = matWADJE(idx2,1);

if isempty(idx3)
    idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else
    idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE
end

% dgamma1 = dgamma2
% B + 2*C*eta = B - 2*C*eta
% !!!!ACROSS A SPLIT, DGAMMA IS CONSTANT FOR ALL PANELS!!!!

% dgamma = B + 2*C*eta
dgamma1 = [ones(len,1) 2.*vecWDVEHVSPN(idx3)];

% dgamma2 = B - 2*C*eta
% Multiplied by -1 as we bring it to the other side of the equal sign
dgamma2 = [ones(len,1) -2.*vecWDVEHVSPN(idx5)].*-1;

% Getting appropriate row and column numbers to assign the above
% dgamma1 and dgamma2 into the D-matrix
col1 = reshape((repmat((idx3.*2)-1,1,2) + repmat(uint32((0:1)), len,1))',[],1);
col2 = reshape((repmat((idx5.*2)-1,1,2) + repmat(uint32((0:1)), len,1))',[],1);
rows = reshape(repmat((1:len)',1,2)',[],1);

% vort = zeros(len,nelements*2);
vort = sparse(len, nelements*2);

vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);

r_vort = zeros(len,1);

% gamma1 = gamma2
% A + B*eta + C*eta^2 = A - B*eta + C*eta^2
% For split:
% gamma12 + gamma22 + ... + gamma2N = gamma41 + gamma42 + ... + gamma4N
% where the left side of the equal side is DVEs with edge 2 at the split
% and the right side is DVEs with edge 4 at the split

% Firstly, the splits are dealt with

% DVE numbers with local edge 2 at the split, along with the degree of the split
[idx13, idx14] = histcounts(idx3, all_DVEs);
dsplit2 = idx14(idx13 > 1);

% Finding the corresponding DVEs with local edge 4 at the split
if ~isempty(dsplit2)
    [idx26,~] = find(repmat(matWADJE(idx2,3),1,length(dsplit2)) == repmat(dsplit2',length(matWADJE(idx2,3)),1));
else
    idx26 = [];
end

% Reassigning dsplit2 so it lines up with dsplit4
dsplit2 = matWADJE(idx2,3);
dsplit2 = double(dsplit2(idx26));
dsplit4 = double(idx4(idx26));

% circ = zeros(len,nelements*2);
circ = sparse(len, nelements*2);
count = 1;

udsplit2 = unique(dsplit2);
for i = 1:length(udsplit2)
    
    idx = find(dsplit2 == udsplit2(i));
    len2 = length(idx);
    % DVE on left of split
    gamma1s = [vecWDVEHVSPN(udsplit2(i)) (2/3).*vecWDVEHVSPN(udsplit2(i)).^2];
    % DVEs on right of split
    gamma2s = [-vecWDVEHVSPN(dsplit4(idx)) (2/3).*vecWDVEHVSPN(dsplit4(idx)).^2].*-1;
    
    col1 = reshape((repmat((udsplit2(i).*2)-1,1,2) + repmat((0:1), 1, 1))',[],1);
    col2 = reshape((repmat((dsplit4(idx).*2)-1,1,2) + repmat((0:1), len2, 1))',[],1);
    
    circ(sub2ind(size(circ),repmat(count,length(col1),1),col1)) = reshape(gamma1s',[],1);
    circ(sub2ind(size(circ),repmat(count,length(col2),1),col2)) = reshape(gamma2s',[],1);
    
    r_circ(count,1) = sum(vecWKGAM(dsplit4(idx))) - vecWKGAM(udsplit2(i));
    
    count = count + 1;
end

% Now we move on to regular circulation between two neighbouring panels
% gamma1 = gamma2
% A + B*eta + C*eta^2 = A - B*eta + C*eta^2

if isempty(idx14(idx13 == 1))
    d2202 = double.empty(0,1);
    d2204 = double.empty(0,1);
else
    % DVE numbers with local edge 2 neighbouring another DVE
    d2202 = idx14(idx13 == 1);
    
    % Finding the DVE numbers with local edge 4 neighbouring the above corresponding DVEs
    d2204 = matWADJE(matWADJE(:,2)==2&ismember(matWADJE(:,1),d2202),3);
end

len3 = length(d2204);

gamma1 = [vecWDVEHVSPN(d2202) (2/3).*vecWDVEHVSPN(d2202).^2];
gamma2 = [-vecWDVEHVSPN(d2204) (2/3).*vecWDVEHVSPN(d2204).^2].*-1;

% Getting appropriate row and column numbers to assign the above
% gamma1 and gamma2 into the D-matrix
col1 = reshape((repmat((d2202.*2)-1,1,2) + repmat((0:1), len3,1))',[],1);
col2 = reshape((repmat((d2204.*2)-1,1,2) + repmat(uint32((0:1)), len3,1))',[],1);
rows = reshape(repmat((1:len3)',1,2)',[],1) + count-1;

circ(sub2ind(size(circ),rows,col1)) = reshape(gamma1',[],1);
circ(sub2ind(size(circ),rows,col2)) = reshape(gamma2',[],1);

r_circ(unique(rows),1) = vecWKGAM(d2204) - vecWKGAM(d2202);

circ(~any(circ,2),:) = [];

% END vorticity and circulation equations between DVEs ----------------------------------------------------------------

% Vorticity equations at symmetry -------------------------------------------------------------------------------------
vort_sym = [];
r_vort_sym = [];
if any(vecWDVESYM) == 1
    
    % DVE number where we have a symmetry condition
    [idx10, ~] = find(vecWDVESYM);
    len = length(idx10);
    locedge = vecWDVESYM(idx10);
    
    dgamma_sym = zeros(len,2);
    
    % dgamma_sym = 0
    % B + 2*C*eta = 0 for right edge (local edge 2), though I doubt this one will be used
    % B - 2*C*eta = 0 for left edge (local edge 4)
    
    dgamma_sym(locedge == 2,:) = [ones(length(idx10(locedge == 2)),1) 2.*vecWDVEHVSPN(idx10(locedge == 2))];
    dgamma_sym(locedge == 4,:) = [ones(length(idx10(locedge == 4)),1) -2.*vecWDVEHVSPN(idx10(locedge == 4))];
    
    % Getting appropriate row and column numbers to assign the above
    % dgamma1 and dgamma2 into the D-matrix
    col3 = reshape((repmat((idx10.*2)-1,1,2) + repmat((0:1), len,1))',[],1);
    rows = reshape(repmat((1:len)',1,2)',[],1);
    
    %     vort_sym = zeros(len,nelements*2);
    vort_sym = sparse(len, nelements*2);
    
    vort_sym(sub2ind(size(vort_sym),rows,col3)) = reshape(dgamma_sym',[],1);
    
    r_vort_sym = zeros(len,1);
    
end
% END vorticity equations at symmetry -----------------------------------------------------------------------------------

% Circulation equations at wingtip --------------------------------------------------------------------------------------

[tip2,~] = find(vecWDVETIP == 2);
[tip4,~] = find(vecWDVETIP == 4);

if sum(vecN) == 1 && ~any(vecWDVESYM); tip2 = [1:length(vecWDVETIP)]'; tip4 = [1:length(vecWDVETIP)]'; end

if isempty(tip2)
    tip2 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(tip4)
    tip4 = double.empty(0,1);
end

gamma1t = [vecWDVEHVSPN(tip2) (2/3).*vecWDVEHVSPN(tip2).^2];
gamma2t = [vecWDVEHVSPN(tip4) -(2/3).*vecWDVEHVSPN(tip4).^2];
gammat = [gamma1t; gamma2t];

col1 = reshape((repmat((tip2.*2)-1,1,2) + repmat((0:1), length(tip2),1))',[],1);
col2 = reshape((repmat((tip4.*2)-1,1,2) + repmat((0:1), length(tip4),1))',[],1);
col = [col1; col2];
rows = reshape(repmat((1:(length(tip2) + length(tip4)))',1,2)',[],1);

% circ_tip = zeros(length(tip2) + length(tip4), nelements*2);
circ_tip = sparse(length(tip2) + length(tip4), nelements*2);

circ_tip(sub2ind(size(circ_tip),rows,col)) = reshape(gammat',[],1);

r_circ_tip = [-vecWKGAM(tip2); vecWKGAM(tip4)];

% END Circulation equations at wingtip ----------------------------------------------------------------------------------

    matWD = [vort; circ; vort_sym; circ_tip];      
    vecWR = [r_vort; r_circ; r_vort_sym; r_circ_tip];

end

