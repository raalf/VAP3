function [ D ] = elementwise3Ddist( A, B )
%ELEMENTWISE3DDIST Summary of this function goes here
%   Detailed explanation goes here
d = A-B;
D = (d(:,:,1).^2+d(:,:,2).^2+d(:,:,3).^2).^0.5;
end

