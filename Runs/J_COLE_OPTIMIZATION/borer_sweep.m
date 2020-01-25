clc
clear

% vinf = [69 73 75 78]
vinf = 78;
[out, ITER, ITEROUTP] = fcnBASELINE_OBJ2(vinf)

% load('borer_paper_data.mat');
% 
% hFig1 = figure(1);
% clf(1);
% 
% plot(overflow_trans(:,1), overflow_trans(:,2), '-ok', 'LineWidth', 1.2);
% hold on
% plot(starccm_trans(:,1), starccm_trans(:,2), '--s', 'Color', [255 12 255]./255, 'LineWidth', 1.2);
% plot(overflow_turb(:,1), overflow_turb(:,2), ':>', 'Color', [0 191 191]./255, 'LineWidth', 1.2);
% plot(fun3d_turb(:,1), fun3d_turb(:,2), '-.d', 'Color', [222 125 0]./255, 'LineWidth', 1.2);
% 
% 
% load('borer_validation_69.mat');
% CL(1) = ITER.CL(end);
% CD(1) = ITER.CD(end);
% 
% load('borer_validation_73.mat');
% CL(2) = ITER.CL(end);
% CD(2) = ITER.CD(end);
% 
% load('borer_validation_75.mat');
% CL(3) = ITER.CL(end);
% CD(3) = ITER.CD(end);
% 
% load('borer_validation_78.mat');
% CL(4) = ITER.CL(end);
% CD(4) = ITER.CD(end);
% 
% plot(CD, CL, '--<r', 'LineWidth', 1.2)
% hold off
% 
% grid minor
% box on
% axis tight
% 
% legend('OVERFLOW (Transitional)', 'STAR-CCM+ (Transitional)', 'OVERFLOW (Turbulent)', 'FUN3D (Turbulent)', 'HOFW (Fixed Wake, Viscous)', 'Location', 'SouthEast')