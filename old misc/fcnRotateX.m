function [result] = fcnRotateX(v1, alpha)
% Facsimile of "rotateX" from vector_algebra.h
% From FreeWake comments, 
% //transforms vector v1 in new co-system that is rotated by alpha
% //around x-axis (RHS!)


% ANGLES IN DEGREES
result(1) = v1(1);
result(2) = v1(2)*cosd(alpha)+v1(3)*sind(alpha);
result(3) = -v1(2)*sind(alpha)+v1(3)*cosd(alpha);

end

