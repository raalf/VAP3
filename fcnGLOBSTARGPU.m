function [xsi] = fcnGLOBSTARGPU(points, roll, pitch, yaw)%, xsi1, xsi2, xsi3, len)
% Transforms a point from local to global reference frame
% Input in RADIANS

% INPUT:
%   points - n x 3 matrix of (x,y,z) points in local coordinates
%   roll - n x 1 vector of roll (nu) angles in radians
%   pitch - n x 1 vector of pitch (eps) angles in radians
%   yaw - n x 1 vector of yaw (psi) angles in radians
% OUTPUT:
%   xsi - n x 3 matrix of "points" in global reference frame

len = length(points(:,1));

cnu = cos(roll);
snu = sin(roll);
ceps = cos(pitch);
seps = sin(pitch);
cpsi = cos(yaw);
spsi = sin(yaw);

% cnu = arrayfun(@cos, gpuArray(roll));
% snu = arrayfun(@sin, gpuArray(roll));
% ceps = arrayfun(@cos, gpuArray(pitch));
% seps = arrayfun(@sin, gpuArray(pitch));
% cpsi = arrayfun(@cos, gpuArray(yaw));
% spsi = arrayfun(@sin, gpuArray(yaw));


xsi = zeros(len,3,'gpuArray');
xsi(:,1) = points(:,1).*(cpsi.*ceps) + points(:,2).*(cpsi.*seps.*snu+spsi.*cnu) + points(:,3).*(-cpsi.*seps.*cnu+spsi.*snu);
xsi(:,2) = points(:,1).*(-spsi.*ceps) + points(:,2).*(-spsi.*seps.*snu+cpsi.*cnu) + points(:,3).*(spsi.*seps.*cnu+cpsi.*snu);
xsi(:,3) = points(:,1).*(seps) + points(:,2).*(-ceps.*snu) + points(:,3).*(ceps.*cnu);

% xsi1 = points1.*(cpsi.*ceps) + points2.*(cpsi.*seps.*snu+spsi.*cnu) + points3.*(-cpsi.*seps.*cnu+spsi.*snu);
% xsi2 = points1.*(-spsi.*ceps) + points2.*(-spsi.*seps.*snu+cpsi.*cnu) + points3.*(spsi.*seps.*cnu+cpsi.*snu);
% xsi3 = points1.*(seps) + points2.*(-ceps.*snu) + points3.*(ceps.*cnu);

end


