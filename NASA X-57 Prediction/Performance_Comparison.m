clc
clear

load('Results/DESIGN_TRANSITION');
load('Results/OVERFLOW_TRANSITION');
load('Results/STARCCM_TRANSITION');
load('Results/OVERFLOW_TURBULENT');
load('Results/FUN3D_TURBULENT');

load('Results/VAP_STEADY_INVISCID_FIXED.mat','vecCL','vecCDI', 'valDELTIME','vecROTORRPM')
vecCLCONV = fcnTIMEAVERAGE(vecCL, vecROTORRPM, valDELTIME);
vecCDICONV = fcnTIMEAVERAGE(vecCDI, vecROTORRPM, valDELTIME);

hFig37 = figure(37);
clf(37);

plot(DESIGN_TRANSITION(:,1), DESIGN_TRANSITION(:,2),'-ok')
hold on
plot(OVERFLOW_TRANSITION(:,1), OVERFLOW_TRANSITION(:,2),'-.^g')
plot(STARCCM_TRANSITION(:,1), STARCCM_TRANSITION(:,2),':sm')
plot(OVERFLOW_TURBULENT(:,1), OVERFLOW_TURBULENT(:,2),'--xb')
plot(FUN3D_TURBULENT(:,1) ,FUN3D_TURBULENT(:,2),'--*r')
plot(vecCDICONV(:), vecCLCONV(:), '-.>k')
hold off

legend('Optimization Result','OVERFLOW (Transition)','STAR-CCM+ (Transition)','OVERFLOW (Turbulent)','FUN3D (Turbulent)','VAP3.0','Location','SouthEast')

grid minor
box on
axis tight

xlabel('C_D','FontSize',15);
ylabel('C_L','FontSize',15);