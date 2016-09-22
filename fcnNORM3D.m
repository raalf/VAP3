function [ M_norm ] = fcnNORM3D( M )
%NORMALIZE3D Summary of this function goes here
%   Detailed explanation goes here


M_length = repmat(((M(:,:,1).^2+M(:,:,2).^2+M(:,:,3).^2).^0.5),1,1,3);
M_norm = M./M_length;



end

