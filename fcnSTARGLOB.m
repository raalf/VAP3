function [xsi] = fcnSTARGLOB(points, roll, pitch, yaw)
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
rot = zeros(3,3,len);

xsi = zeros(len,1,3);

x = permute(points,[1 3 2]);

[A,B,~] = size(x);

cnu = cos(roll);
snu = sin(roll);
ceps = cos(pitch);
seps = sin(pitch);
cpsi = cos(yaw);
spsi = sin(yaw);

% rot(1,1,:) = cpsi.*ceps;
% rot(1,2,:) = cpsi.*seps.*snu+spsi.*cnu;
% rot(1,3,:) = -cpsi.*seps.*cnu+spsi.*snu;
% rot(2,1,:) = -spsi.*ceps;
% rot(2,2,:) = -spsi.*seps.*snu+cpsi.*cnu;
% rot(2,3,:) = spsi.*seps.*cnu+cpsi.*snu;
% rot(3,1,:) = seps;
% rot(3,2,:) = -ceps.*snu;
% rot(3,3,:) = ceps.*cnu;

% xsi(:,:,1) = x(:,:,1).*reshape(rot(1,1,:),A,B) + x(:,:,2).*reshape(rot(2,1,:),A,B) + x(:,:,3).*reshape(rot(3,1,:),A,B);
% xsi(:,:,2) = x(:,:,1).*reshape(rot(1,2,:),A,B) + x(:,:,2).*reshape(rot(2,2,:),A,B) + x(:,:,3).*reshape(rot(3,2,:),A,B);
% xsi(:,:,3) = x(:,:,1).*reshape(rot(1,3,:),A,B) + x(:,:,2).*reshape(rot(2,3,:),A,B) + x(:,:,3).*reshape(rot(3,3,:),A,B);

xsi(:,:,1) = x(:,:,1).*reshape(cpsi.*ceps,A,B) + x(:,:,2).*reshape(-spsi.*ceps,A,B) + x(:,:,3).*reshape(seps,A,B);
xsi(:,:,2) = x(:,:,1).*reshape(cpsi.*seps.*snu+spsi.*cnu,A,B) + x(:,:,2).*reshape(-spsi.*seps.*snu+cpsi.*cnu,A,B) + x(:,:,3).*reshape(-ceps.*snu,A,B);
xsi(:,:,3) = x(:,:,1).*reshape(-cpsi.*seps.*cnu+spsi.*snu,A,B) + x(:,:,2).*reshape(spsi.*seps.*cnu+cpsi.*snu,A,B) + x(:,:,3).*reshape(ceps.*cnu,A,B);

xsi = reshape(xsi,[],3,1);

end

