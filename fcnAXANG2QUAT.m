function q = fcnAXANG2QUAT(axang)

% For a single axis-angle vector [ax ay az t] the output quaternion
% q = [w x y z], can be computed as follows:
% w = cos(t/2)
% x = ax*sin(t/2)
% y = ay*sin(t/2)
% z = az*sin(t/2)

v = bsxfun(@times, axang(:,1:3), 1./sqrt(sum(axang(:,1:3).^2,2)));

% Create the quaternion
thetaHalf = axang(:,4)/2;
sinThetaHalf = sin(thetaHalf);
q = [cos(thetaHalf), v(:,1).*sinThetaHalf, v(:,2).*sinThetaHalf, v(:,3).*sinThetaHalf];

end

