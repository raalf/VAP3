function [matE] = fcnWEWAKE(valNELE, matADJE, vecDVEHVCRD)
% OPERA LITE - CHORDWISE CIRCULATION AND VORTICITY

all_DVEs = [1:valNELE]';

% POSITIVE XSI DIRECTION IS REARWARD, ALIGNED WITH GLOBAL X

% Vorticity and circulation equations between DVEs ----------------------------------------------------------------------

idx1 = matADJE(:,2) == 1; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 3; % Borders on 4th local edge (left side)

len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
% This returns a vector where we can see for edge 2 equations, which is the
% matching edge 4 equation
% This will work well for split vorticity but it will need modification for split circulation later on
[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);

if isempty(idx3)
    idx5 = idx3; % Ensuring the empty is the correct size if this is empty
else
    idx5 = idx4(eqn_num); % Corresponding edge 3 DVE for every edge 1 DVE
end

% dgamma1 = dgamma2
% B + 2*C*eta = B - 2*C*eta
% !!!!ACROSS A SPLIT, DGAMMA IS CONSTANT FOR ALL PANELS!!!!

cols_zeros = zeros(len,2);
% dgamma = B + 2*C*eta
dgamma1 = [cols_zeros ones(len,1) -2.*vecDVEHVCRD(idx3)];

% dgamma2 = B - 2*C*eta
% Multiplied by -1 as we bring it to the other side of the equal sign
dgamma2 = [cols_zeros ones(len,1) 2.*vecDVEHVCRD(idx5)].*-1;

% Getting appropriate row and column numbers to assign the above
% dgamma1 and dgamma2 into the D-matrix
col1 = reshape([repmat((idx3.*4)-3,1,4) + repmat([0:3], len,1)]',[],1);
col2 = reshape([repmat((idx5.*4)-3,1,4) + repmat([0:3], len,1)]',[],1);
rows = reshape(repmat([1:len]',1,4)',[],1);

% vort = zeros(len,valNELE*3);
vort = sparse(len,valNELE*4);

vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);


% [idx13, idx14] = histcounts(idx3, all_DVEs);
% 
% % circ = zeros(len,valNELE*3);
% circ = sparse(len,valNELE*4);
% 
% % Now we move on to regular circulation between two neighbouring panels
% % gamma1 = gamma2
% % A + B*eta + C*eta^2 = A - B*eta + C*eta^2
% 
% if isempty(idx14(idx13 == 1))
%     d2202 = double.empty(0,1);
%     d2204 = double.empty(0,1);
% else
%     % DVE numbers with local edge 2 neighbouring another DVE
%     d2202 = idx14(idx13 == 1);
%     
%     % Finding the DVE numbers with local edge 3 neighbouring the above corresponding DVEs
%     [idx27,~] = find(repmat(matADJE(idx1,1),1,length(d2202)) == repmat(d2202',length(matADJE(idx1,1)),1));
%     d2204 = matADJE(matADJE(:,2) == 1,3);
%     d2204 = d2204(idx27);
% end
% 
% 
% %%%%% WHY DO I NEED THE ABOVE D2202 AND D2204 STUFF??? SPLIT WING?
% d2202 = idx3;
% d2204 = idx4; 
% 
% len3 = length(d2204);
% 
% gamma1 = [zeros(len3,2) -vecDVEHVCRD(d2202) vecDVEHVCRD(d2202).^2];
% gamma2 = [zeros(len3,2) vecDVEHVCRD(d2204) vecDVEHVCRD(d2204).^2].*-1;
% 
% % Getting appropriate row and column numbers to assign the above
% % gamma1 and gamma2 into the D-matrix
% col1 = reshape([repmat((d2202.*4)-3,1,4) + repmat([0:3], len3,1)]',[],1);
% col2 = reshape([repmat((d2204.*4)-3,1,4) + repmat([0:3], len3,1)]',[],1);
% rows = reshape(repmat([1:len3]',1,4)',[],1);
% 
% circ(sub2ind(size(circ),rows,col1)) = reshape(gamma1',[],1);
% circ(sub2ind(size(circ),rows,col2)) = reshape(gamma2',[],1);
% 
% circ(~any(circ,2),:) = [];

% END vorticity and circulation equations between DVEs ----------------------------------------------------------------

% Vorticity equations at leading and trailing edges -----------------------------------------------------------------

% [edge1,~] = find(vecDVELE == 1 | vecDVETE == 1);
% [edge3,~] = find(vecDVELE == 3 | vecDVETE == 3);

[edge1,~] = find(vecDVETE == 1);
[edge3,~] = find(vecDVETE == 3);

if isempty(edge1)
    edge1 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
elseif isempty(edge3)
    edge3 = double.empty(0,1);
end

dgamma1t = [ones(length(edge1),1) -2.*vecDVEHVCRD(edge1)];
dgamma2t = [ones(length(edge3),1) 2.*vecDVEHVCRD(edge3)];
% dgamma1t = [ones(length(edge1),1) zeros(length(edge1),2) vecDVEHVCRD(edge1) vecDVEHVCRD(edge1).^2];
% dgamma2t = [ones(length(edge3),1) zeros(length(edge3),2) -vecDVEHVCRD(edge3) vecDVEHVCRD(edge3).^2];
dgammat = [dgamma1t; dgamma2t];

col1 = reshape([repmat((edge1.*4)-3,1,4) + repmat([0:3], length(edge1),1)]',[],1);
col2 = reshape([repmat((edge3.*4)-3,1,4) + repmat([0:3], length(edge3),1)]',[],1);
col = [col1; col2];
rows = reshape(repmat([1:(length(edge1) + length(edge3))]',1,4)',[],1);

% circ_tip = zeros(length(tip2) + length(tip4), valNELE*3);
circ_freeedge = sparse(length(edge1) + length(edge3), valNELE*4);

circ_freeedge(sub2ind(size(circ_freeedge),rows,col)) = reshape(dgammat',[],1);

% END Circulation equations at wingtip ----------------------------------------------------------------------------------

matE = [vort; circ; circ_freeedge];

end

