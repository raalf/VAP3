clc
clear

%%
load('Results/DESIGN_TRANSITION');
load('Results/OVERFLOW_TRANSITION');
load('Results/STARCCM_TRANSITION');
load('Results/OVERFLOW_TURBULENT');
load('Results/FUN3D_TURBULENT');

load('Results/VAP_STEADY_VISCOUS_FIXED.mat','vecCLv','vecCD', 'valDELTIME','vecROTORRPM')
% vecCLCONV = fcnTIMEAVERAGE(vecCL, vecROTORRPM, valDELTIME);
% vecCDICONV = fcnTIMEAVERAGE(vecCDI, vecROTORRPM, valDELTIME);

hFig37 = figure(37);
clf(37);

% plot(DESIGN_TRANSITION(:,1), DESIGN_TRANSITION(:,2),'-ok')
hold on
plot(OVERFLOW_TRANSITION(:,1), OVERFLOW_TRANSITION(:,2),'-.^g')
plot(STARCCM_TRANSITION(:,1), STARCCM_TRANSITION(:,2),':sm')
% plot(OVERFLOW_TURBULENT(:,1), OVERFLOW_TURBULENT(:,2),'--xb')
% plot(FUN3D_TURBULENT(:,1) ,FUN3D_TURBULENT(:,2),'--*r')
plot(vecCD, vecCLv(:), '-.>k')
hold off

% legend('Optimization Objective Function','OVERFLOW (Transitional)','STAR-CCM+ (Transitional)','OVERFLOW (Turbulent)','FUN3D (Turbulent)','VAP3.0','Location','SouthEast')
legend('OVERFLOW (Transitional)','STAR-CCM+ (Transitional)','VAP3.0','Location','SouthEast')

grid minor
box on
axis tight

xlabel('C_D','FontSize',15);
ylabel('C_L','FontSize',15);

%% Propeller Prediction
clear

load('Results/NASA_X57_CRUISE_PROP_BvT.mat')
load('Results/VAP_NASA_X57_CRUISE_PROP.mat','vecCT','vecROTORRPM','vecROTDIAM','vecCOLLECTIVE')

thrust = vecCT(end,:,:).*(((vecROTORRPM/60).^2).*((vecROTDIAM).^4));
thrust = thrust(:);

hFig38 = figure(38);
clf(38);

plot(NASA_X57_CRUISE_PROP_BvT(:,1),NASA_X57_CRUISE_PROP_BvT(:,2),'-ok');
hold on
plot(vecCOLLECTIVE + 67, thrust, '--sb') % 67 is our root twist
hold off

legend('XROTOR (Patterson)','VAP3.0 (Inviscid)','Location','NorthWest')

grid minor
box on
axis tight

xlabel('Blade Root (\beta) Angle (Degrees)','FontSize',15);
ylabel('Thrust (N)','FontSize',15);