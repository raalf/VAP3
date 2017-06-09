function [xo, u, normal] = fcnNewxo(FW, x_left, x_right, xo_old)
% 	//This function computes a new reference point and the vector
% 	//between left and right edge
% 	//
% 	//input:
% 	//	info		- general information
% 	//	x_left		- point to the left
% 	//	x_right		- point to the right
% 	//	wakeDVE		- information of wake DVE of interest
% 	//
% 	//ouput:
% 	// 	For the wake DVE of interest the following values are being updated
% 	//	x[]			- reference point in global ref.frame
% 	//	u[]			- velocity in x[]

normal = 0.5*(x_right-x_left); % Storing this temporarily (it is vector from left to right side)

xo = x_left + normal;

u = (xo - xo_old)/FW.Deltime;

end

