function [matWE] = fcnWEWAKE(all_DVEs, matWADJE, vecWDVEHVCRD, vecWDVELE, vecWDVETE, matCOEFF, vecDVEHVCRD, vecDVETE, vecWEKGAM)
% OPERA LITE - CHORDWISE CIRCULATION AND VORTICITY

nelements = length(all_DVEs);

% POSITIVE XSI DIRECTION IS REARWARD, ALIGNED WITH GLOBAL X

% Vorticity and circulation equations between DVEs ----------------------------------------------------------------------

idx1 = matWADJE(:,2) == 1; % Borders on 2nd local edge (right side)
idx2 = matWADJE(:,2) == 3; % Borders on 4th local edge (left side)

len = length(idx1(idx1>0,1));

if len > 0
    
    % Mapping equations evaluated on edge 2 to equations evaluated on edge 4
    % This returns a vector where we can see for edge 2 equations, which is the
    % matching edge 4 equation
    % This will work well for split vorticity but it will need modification for split circulation later on
    [~, LB] = ismembertol([matWADJE(idx1,1) matWADJE(idx1,3)],[matWADJE(idx2,3) matWADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
    eqn_num = cell2mat(LB);
    
    idx3 = matWADJE(idx1,1);
    idx4 = matWADJE(idx2,1);
    
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
    dgamma1 = [cols_zeros ones(len,1) -2.*vecWDVEHVCRD(idx3)];
    
    % dgamma2 = B - 2*C*eta
    % Multiplied by -1 as we bring it to the other side of the equal sign
    dgamma2 = [cols_zeros ones(len,1) 2.*vecWDVEHVCRD(idx5)].*-1;
    
    % Getting appropriate row and column numbers to assign the above
    % dgamma1 and dgamma2 into the D-matrix
    col1 = reshape([repmat((idx3.*4)-3,1,4) + repmat([0:3], len,1)]',[],1);
    col2 = reshape([repmat((idx5.*4)-3,1,4) + repmat([0:3], len,1)]',[],1);
    rows = reshape(repmat([1:len]',1,4)',[],1);
    
    % vort = zeros(len,valNELE*3);
    vort = sparse(len,nelements*4);
    
    vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
    vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);
    
    r_vort = zeros(len,1);
    
    
    [idx13, idx14] = histcounts(idx3, all_DVEs);
    
    % circ = zeros(len,valNELE*3);
    circ = sparse(len,nelements*4);
    
    % Now we move on to regular circulation between two neighbouring panels
    % gamma1 = gamma2
    % A + B*eta + C*eta^2 = A - B*eta + C*eta^2
    
    if isempty(idx14(idx13 == 1))
        d2202 = double.empty(0,1);
        d2204 = double.empty(0,1);
    else
        % DVE numbers with local edge 2 neighbouring another DVE
        d2202 = idx14(idx13 == 1);
        
        % Finding the DVE numbers with local edge 3 neighbouring the above corresponding DVEs
        [idx27,~] = find(repmat(matWADJE(idx1,1),1,length(d2202)) == repmat(d2202',length(matWADJE(idx1,1)),1));
        d2204 = matWADJE(matWADJE(:,2) == 1,3);
        d2204 = d2204(idx27);
    end
    
    
    %%%%% WHY DO I NEED THE ABOVE D2202 AND D2204 STUFF??? SPLIT WING?
    d2202 = idx3;
    d2204 = idx4;
    
    len3 = length(d2204);
    
    gamma1 = [cols_zeros -vecWDVEHVCRD(d2202) (2/3).*vecWDVEHVCRD(d2202).^2];
    gamma2 = [cols_zeros vecWDVEHVCRD(d2204) (2/3).*vecWDVEHVCRD(d2204).^2].*-1;
    
    % Getting appropriate row and column numbers to assign the above
    % gamma1 and gamma2 into the D-matrix
    col1 = reshape([repmat((d2202.*4)-3,1,4) + repmat([0:3], len3,1)]',[],1);
    col2 = reshape([repmat((d2204.*4)-3,1,4) + repmat([0:3], len3,1)]',[],1);
    rows = reshape(repmat([1:len3]',1,4)',[],1);
    
    circ(sub2ind(size(circ),rows,col1)) = reshape(gamma1',[],1);
    circ(sub2ind(size(circ),rows,col2)) = reshape(gamma2',[],1);
    
    r_circ = zeros(size(gamma1,1),1);
    r_circ(unique(rows),1) = vecWEKGAM(d2204) - vecWEKGAM(d2202); % RIGHT OR WRONG WAY ROUND?
    
    circ(~any(circ,2),:) = [];
   
    
else
    vort = [];
    r_vort = [];
    circ = [];
    r_circ = [];
end

% END vorticity and circulation equations between DVEs ----------------------------------------------------------------

% Vorticity equations at trailing edge of WAKE -----------------------------------------------------------------

[edge1,~] = find(vecWDVETE > 0);

if isempty(edge1)
    edge1 = double.empty(0,1); % Ensuring the empty is the correct size if this is empty
end

dgamma1t = [zeros(length(edge1),2) ones(length(edge1),1) 2.*vecWDVEHVCRD(edge1)];
% dgamma1t = [ones(length(edge1),1) zeros(length(edge1),2) vecDVEHVCRD(edge1) vecDVEHVCRD(edge1).^2];


col = reshape([repmat((edge1.*4)-3,1,4) + repmat([0:3], length(edge1),1)]',[],1);
rows = reshape(repmat([1:length(edge1)]',1,4)',[],1);

vort_te = sparse(length(edge1), nelements*4);

vort_te(sub2ind(size(vort_te),rows,col)) = reshape(dgamma1t,[],1);

r_vort_te = zeros(length(edge1),1);


% Circulation at leading edge of WAKE --------------------------------------------------------

[edge3,~] = find(vecWDVELE == 1);

if isempty(edge3)
    edge3 = double.empty(0,1);
end

gamma2t = [zeros(length(edge3),2) ones(length(edge3),1) -(2/3).*vecWDVEHVCRD(edge3).^2];

col = reshape([repmat((edge3.*4)-3,1,4) + repmat([0:3], length(edge3),1)]',[],1);
rows = reshape(repmat([1:length(edge3)]',1,4)',[],1);

% circ_tip = zeros(length(tip2) + length(tip4), valNELE*3);
circ_wle = sparse(length(edge3), nelements*4);

circ_wle(sub2ind(size(circ_wle),rows,col)) = reshape(gamma2t,[],1);

r_circ_le = zeros(length(edge3),1);

r_circ_le = 1.*[matCOEFF(vecDVETE > 0,4) + (2/3).*matCOEFF(vecDVETE > 0,5).*vecDVEHVCRD(vecDVETE > 0).^2];


matWE = [vort; circ; vort_te; circ_wle];
vecRE = [r_vort; r_circ; r_vort_te; r_circ_le];


end

