function [matD] = fcnDEXPAND_OL(matD, matE, valNELE)
% This function expands the D-matrix or WD-matrix
% to combine both D and E matrices in OPERA-lite
% valNELE = valNELE for surface, valWNELE for wake

cols = size(matD,2)/valNELE; % Should be 3 for surface, 2 for wake

matDxpand = zeros(valNELE*2, valNELE*(2+cols));

% Getting old column numbers and new column numbers
old_col = [1:valNELE*cols];
new_col = [1:2:valNELE*2]-1;
new_col = old_col + reshape(repmat(new_col',1,cols)',[],1)';

matDxpand(:,new_col) = matD;

matD = [matDxpand; matE];


end

