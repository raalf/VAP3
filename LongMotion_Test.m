close all
dt = 0.001;
Tmax = 50;
TRIM.matA = [-0.0008 1.8365 -0.4005 -9.81; -0.685 -10.7835 -0.5031 0;...
    0.03 -0.4574 -17.3427 0; 0 0 1 0];
[TRIM.vecEIGVEC, TRIM.vecEIGVAL] = eig(TRIM.matA);
TRIM.vecEIGVAL = diag(TRIM.vecEIGVAL);

i = 1;

for t = 0:dt:Tmax
    
    for j = 1:size(TRIM.matA,1)
        x_t(:,j) = TRIM.vecEIGVEC(:,j)*exp(TRIM.vecEIGVAL(j).*t);
    end

    xt = sum(x_t,2);

    delu(i) = xt(1);
    w(i) = xt(2);
    q(i) = xt(3);
    deltheta(i) = xt(4);
    
    i = i + 1;
    
end

% [x0, sys_eig] = eig(XFLR_matA); % Find eigenvectors and eigenvalues of the system
% sys_eig = diag(sys_eig);
% 
% i = 1;
% for t = 0:dt:Tmax
%     
%     for j = 1:size(XFLR_matA,1)
%         x_t(:,j) = x0(:,j)*exp(sys_eig(j).*t);
%     end
% 
%     xt = sum(x_t,2);
% 
%     delu_xflr(i) = xt(1);
%     w_xflr(i) = xt(2);
%     q_xflr(i) = xt(3);
%     deltheta_xflr(i) = xt(4);
%     
%     i = i + 1;
%     
% end


figure(69)
plot(0:dt:Tmax,delu,'-b','linewidth',1.5)
hold on
plot(0:dt:Tmax,w,'-r','linewidth',1.5)
% plot(0:dt:Tmax,delu_xflr,'--b','linewidth',1.5)
% plot(0:dt:Tmax,w_xflr,'--r','linewidth',1.5)
grid on
xlabel('Time (s)')
ylabel('Perturbation (m/s)')
legend('{\Delta}u','w','location','southeast')

figure(420)
plot(0:dt:Tmax,q,'-b','linewidth',1.5)
hold on
plot(0:dt:Tmax,deltheta,'-r','linewidth',1.5)
grid on
% plot(0:dt:Tmax,q_xflr,'--b','linewidth',1.5)
% plot(0:dt:Tmax,deltheta_xflr,'--r','linewidth',1.5)
xlabel('Time (s)')
ylabel('Perturbation (deg, deg/s)')
legend('q','{\Delta}{\theta}','location','northeast')