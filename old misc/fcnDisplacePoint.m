function [x] = fcnDisplacePoint(FW, uk, ul, um, x)
% Copy of Displace_Wake from wake_geometry.cpp

% 	//	function displaces a point along the local streamline
% 	//
% 	//input:
% 	//	info		- general information
% 	//	rightnow	- current time step
% 	//	uk			- local velocity at midchord of DVE downstream
% 	//	ul			- local velocity at midchord of DVE of interest
% 	//	um			- local velocity at midchord of DVE upstream
% 	//  x			- point that is displaced
% 	//
% 	//ouput:
% 	// 	For each wake DVE in wakePtr updated
% 	//	x[]			- point of interest in global ref.frame

u = 0.25.*(uk + 2.*ul + um);

x = x + (u.*FW.Deltime*FW.Uinf);


end

