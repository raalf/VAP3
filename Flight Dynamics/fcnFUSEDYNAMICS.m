function dydt = fcnFUSEDYNAMICS(t,y,X,Z,M,mtot,mf,g,Iyy)

dydt = zeros(4,1);

dydt(1) = X/mtot - g*sin(y(4)) - y(3)*y(2);
dydt(2) = Z/mf + g*cos(y(4)) + y(3)*y(1);
dydt(3) = M/Iyy;
dydt(4) = y(3);
