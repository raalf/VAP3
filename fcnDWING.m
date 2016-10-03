function [matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP)
% Currently creates the upper 2/3rds of the D-matrix, which are the boundary conditions based
% on adjacent DVEs, etc.

% INPUT:
%   valNELE - total number of DVEs
%   matADJE - ? x 3 adjacency matrix, where columns are: DVE | local edge | adjacent DVE
%   vecDVEHVSPN - valNELE x 1 vector of DVE half spans
%   vecDVESYM - valNELE x 1 vector of which DVEs have symmetry on which edge (0 for no symmetry,
%                   2 for local edge 2, 4 for local edge 4)
%   vecDVETIP - valNELE x 1 vector of which DVEs are at the wingtip. Similar format to vecDVESYM above
% OUTPUT:
%   matD - currently, upper 2/3 of D-matrix

all_DVEs = [1:valNELE]';

% Vorticity and circulation equations between DVEs ----------------------------------------------------------------------

idx1 = matADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 4; % Borders on 4th local edge (left side)

len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
% This returns a vector where we can see for edge 2 equations, which is the
% matching edge 4 equation
% This will work well for split vorticity but it will need modification for split circulation later on
% [~, LB] = ismembertol([matADJE(idx2,3) matADJE(idx2,1)],[matADJE(idx1,1) matADJE(idx1,3)],'ByRows',true,'OutputAllIndices',true);
[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);
idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE

% dgamma1 = dgamma2
% B + 2*C*eta = B - 2*C*eta
% !!!!ACROSS A SPLIT, DGAMMA IS CONSTANT FOR ALL PANELS!!!!

% dgamma = B + 2*C*eta
dgamma1 = [zeros(len,1) ones(len,1) 2.*vecDVEHVSPN(idx3)];

% dgamma2 = B - 2*C*eta
% Multiplied by -1 as we bring it to the other side of the equal sign
dgamma2 = [zeros(len,1) ones(len,1) -2.*vecDVEHVSPN(idx5)].*-1;

% Getting appropriate row and column numbers to assign the above
% dgamma1 and dgamma2 into the D-matrix
col1 = reshape([repmat((idx3.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
col2 = reshape([repmat((idx5.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
rows = reshape(repmat([1:len]',1,3)',[],1);

% vort = zeros(len,valNELE*3);
vort = sparse(len,valNELE*3);

vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);

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
    [idx26,~] = find(repmat(matADJE(idx2,3),1,length(dsplit2)) == repmat(dsplit2',length(matADJE(idx2,3)),1));
else
    idx26 = [];
end
% Reassigning dsplit2 so it lines up with dsplit4
dsplit2 = matADJE(idx2,3);
dsplit2 = dsplit2(idx26);
dsplit4 = idx4(idx26);

% circ = zeros(len,valNELE*3);
circ = sparse(len,valNELE*3);

count = 1;

udsplit2 = unique(dsplit2);
for i = 1:length(udsplit2)
    
    idx = find(dsplit2 == udsplit2(i));
    len2 = length(idx);
    % DVE on left of split
    gamma1s = [1 vecDVEHVSPN(udsplit2(i)) vecDVEHVSPN(udsplit2(i)).^2];
    % DVEs on right of split
    gamma2s = [ones(len2,1) -vecDVEHVSPN(dsplit4(idx)) vecDVEHVSPN(dsplit4(idx)).^2].*-1;
    
    col1 = reshape([repmat((udsplit2(i).*3)-2,1,3) + repmat([0:2], 1, 1)]',[],1);
    col2 = reshape([repmat((dsplit4(idx).*3)-2,1,3) + repmat([0:2], len2, 1)]',[],1);
    
    circ(sub2ind(size(circ),repmat(count,length(col1),1),col1)) = reshape(gamma1s',[],1);
    circ(sub2ind(size(circ),repmat(count,length(col2),1),col2)) = reshape(gamma2s',[],1);
    
    count = count + 1;
end

% Now we move on to regular circulation between two neighbouring panels
% gamma1 = gamma2
% A + B*eta + C*eta^2 = A - B*eta + C*eta^2

% DVE numbers with local edge 2 neighbouring another DVE
d2202 = idx14(idx13 == 1);

% Finding the DVE numbers with local edge 4 neighbouring the above corresponding DVEs
[idx27,~] = find(repmat(matADJE(idx1,1),1,length(d2202)) == repmat(d2202',length(matADJE(idx1,1)),1));
d2204 = matADJE(matADJE(:,2) == 2,3);
d2204 = d2204(idx27);

len3 = length(d2204);

gamma1 = [ones(len3,1) vecDVEHVSPN(d2202) vecDVEHVSPN(d2202).^2];
gamma2 = [ones(len3,1) -vecDVEHVSPN(d2204) vecDVEHVSPN(d2204).^2].*-1;

% Getting appropriate row and column numbers to assign the above
% gamma1 and gamma2 into the D-matrix
col1 = reshape([repmat((d2202.*3)-2,1,3) + repmat([0:2], len3,1)]',[],1);
col2 = reshape([repmat((d2204.*3)-2,1,3) + repmat([0:2], len3,1)]',[],1);
rows = reshape(repmat([1:len3]',1,3)',[],1) + count-1;

circ(sub2ind(size(circ),rows,col1)) = reshape(gamma1',[],1);
circ(sub2ind(size(circ),rows,col2)) = reshape(gamma2',[],1);

circ(~any(circ,2),:) = [];

% END vorticity and circulation equations between DVEs ----------------------------------------------------------------

% Vorticity equations at symmetry -------------------------------------------------------------------------------------
if ~isempty(vecDVESYM) == 1
    
    % DVE number where we have a symmetry condition
    [idx10, ~] = find(vecDVESYM);
    len = length(idx10);
    locedge = vecDVESYM(idx10);
    
    dgamma_sym = zeros(len,3);
    
    % dgamma_sym = 0
    % B + 2*C*eta = 0 for right edge (local edge 2), though I doubt this one will be used
    % B - 2*C*eta = 0 for left edge (local edge 4)
    
    dgamma_sym(locedge == 2,:) = [zeros(length(idx10(locedge == 2)),1) ones(length(idx10(locedge == 2)),1) 2.*vecDVEHVSPN(idx10(locedge == 2))];
    dgamma_sym(locedge == 4,:) = [zeros(length(idx10(locedge == 4)),1) ones(length(idx10(locedge == 4)),1) -2.*vecDVEHVSPN(idx10(locedge == 4))];
    
    % Getting appropriate row and column numbers to assign the above
    % dgamma1 and dgamma2 into the D-matrix
    col3 = reshape([repmat((idx10.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
    rows = reshape(repmat([1:len]',1,3)',[],1);
    
    %     vort_sym = zeros(len,valNELE*3);
    vort_sym = sparse(len,valNELE*3);
    
    vort_sym(sub2ind(size(vort_sym),rows,col3)) = reshape(dgamma_sym',[],1);
    
end
% END vorticity equations at symmetry -----------------------------------------------------------------------------------

% Circulation equations at wingtip --------------------------------------------------------------------------------------

[tip2,~] = find(vecDVETIP == 2);
[tip4,~] = find(vecDVETIP == 4);

gamma1t = [ones(length(tip2),1) vecDVEHVSPN(tip2) vecDVEHVSPN(tip2).^2];
gamma2t = [ones(length(tip4),1) -vecDVEHVSPN(tip4) vecDVEHVSPN(tip4).^2];
gammat = [gamma1t; gamma2t];

col1 = reshape([repmat((tip2.*3)-2,1,3) + repmat([0:2], length(tip2),1)]',[],1);
col2 = reshape([repmat((tip4.*3)-2,1,3) + repmat([0:2], length(tip4),1)]',[],1);
col = [col1; col2];
rows = reshape(repmat([1:(length(tip2) + length(tip4))]',1,3)',[],1);

% circ_tip = zeros(length(tip2) + length(tip4), valNELE*3);
circ_tip = sparse(length(tip2) + length(tip4), valNELE*3);

circ_tip(sub2ind(size(circ_tip),rows,col)) = reshape(gammat',[],1);

% END Circulation equations at wingtip ----------------------------------------------------------------------------------

matD = [vort; circ; vort_sym; circ_tip];

end

