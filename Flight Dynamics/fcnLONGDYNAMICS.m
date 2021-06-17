function dydt = fcnLONGDYNAMICS(t,y,X,Z,M,m,mf,g,Iyy)

% System of differential equations for longitudinal motion

% X - Resultant aerodynamic force in the X-direction (positive is forward)
% Z - Resultant aerodynamic force in the Z-direction (positive is down)
% M - Resultant pitching moment (positive is nose up)
% Iyy - Mass moment of inertia about the Y-axis

dydt = zeros(4,1);

dydt(1) = X/m - g*sin(y(4)) - y(3)*y(2); % Acceleration in u
dydt(2) = Z/mf + g*cos(y(4)) + y(3)*y(1); % Acceleration in w
dydt(3) = M/Iyy; % Pitch rate
dydt(4) = y(3); % Pitch angle

end