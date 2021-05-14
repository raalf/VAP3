load('Sharp_Edge_m1.mat')

figure(6969)
hold on
plot((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME,OUTP.vecCL(COND.valGUSTSTART:end),'-','linewidth',1.5)
box on
grid on
ylabel('C_l')
xlabel('Time (s)')

load('Sharp_Edge_m5.mat')

figure(6969)
hold on
plot((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME,OUTP.vecCL(COND.valGUSTSTART:end),'-','linewidth',1.5)
box on
grid on
ylabel('C_l')
xlabel('Time (s)')

load('Sharp_Edge_m10.mat')

figure(6969)
hold on
plot((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME,OUTP.vecCL(COND.valGUSTSTART:end),'-','linewidth',1.5)
box on
grid on
ylabel('C_l')
xlabel('Time (s)')

load('Sharp_Edge_m30.mat')

figure(6969)
hold on
plot((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME,OUTP.vecCL(COND.valGUSTSTART:end),'-','linewidth',1.5)
box on
grid on
ylabel('C_l')
xlabel('Time (s)')

Kussner_Function

plot(t,Cl,'-k','linewidth',1.5)

legend('m = 1','m = 5','m = 10','m = 30','Kussner','location','southeast')