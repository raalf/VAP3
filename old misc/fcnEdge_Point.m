function x = fcnEdge_Point(FW,i,j,phi,eta,xsi,xo)
% phi needs to be in deg
% //edge location in local ref. frame
edge(1) = xsi+eta*tand(phi);
edge(2)=eta;
edge(3)=0;

% //transformation of xsi-edge into global reference frame
x = star_glob( edge,FW.Panels(i).DVE.roll(j),FW.Panels(i).DVE.pitch(j),FW.Panels(i).DVE.yaw(j) );

% //with respect to reference point
x = x+xo;
