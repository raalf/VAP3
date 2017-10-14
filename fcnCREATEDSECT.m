function [sect] = fcnCREATEDSECT(sect, len, vars, first_idx, second_idx, first_part, second_part)
% Getting appropriate row and column numbers to assign the above
% parts into the D-matrix
% vars is number of variables, 3 for VAP and 5 for HDVE

rows = reshape(repmat([1:len]',1,vars)',[],1);
col1 = reshape([repmat((first_idx.*vars)-(vars-1),1,vars) + repmat([0:(vars-1)], len,1)]',[],1);
sect(sub2ind(size(sect),rows,col1)) = reshape(first_part',[],1);

if ~isempty(second_idx)
    col2 = reshape([repmat((second_idx.*vars)-(vars-1),1,vars) + repmat([0:(vars-1)], len,1)]',[],1);
    sect(sub2ind(size(sect),rows,col2)) = reshape(second_part',[],1);
end
end