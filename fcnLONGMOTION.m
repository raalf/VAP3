function [delu, w, q, deltheta] = fcnLONGMOTION(Xu, Xw, Zu, Zw, Zq, Zwdot, Mu, Mw, Mq, Mwdot, m, g, Iy, theta0, u0, t)

% Function to solve longitudinal equations of motion via small disturbance
% theory

A = [Xu/m, Xw/m, 0, -g*cos(theta0);...
    Zu/(m-Zwdot), Zw/(m-Zwdot), (Zq + m*u0)/(m-Zwdot), (-m*g*sin(theta0))/(m-Zwdot);...
    (1/Iy)*(Mu + (Mwdot*Zu)/(m-Zwdot)), (1/Iy)*(Mw + (Mwdot*Zw)/(m-Zwdot)), (1/Iy)*(Mq + (Mwdot*(Zq + m*u0))/(m-Zwdot)), (-Mwdot*m*g*sin(theta0))/(Iy*(m-Zwdot));...
    0, 0, 1, 0];

[x0, sys_eig] = eig(A); % Find eigenvectors and eigenvalues of the system
sys_eig = diag(sys_eig);

for i = 1:size(A,1)
    x_t(:,i) = x0(:,i)*exp(sys_eig(i).*t);
end

xt = sum(x_t,2);

delu = xt(1);
w = xt(2);
q = xt(3);
deltheta = xt(4);
