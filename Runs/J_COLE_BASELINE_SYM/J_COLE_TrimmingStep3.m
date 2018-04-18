clear
clc

% Import 1st Iteration Data
load('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps_2ndTrimIter.mat')
OUTP2 = OUTP;
clear OUTP;


CL2 = nanmean([OUTP2.vecCL],1);
CT2 = nanmean([OUTP2.vecCT],1);

CD1temp = reshape([OUTP1.vecCD],[],length(OUTP1));
CD1temp(isnan([OUTP1.vecCL])) = nan;
CD1 = nanmean(CD1temp,1);

CD2temp = reshape([OUTP2.vecCD],[],length(OUTP2));
CD2temp(isnan([OUTP2.vecCL])) = nan;
CD2 = nanmean(CD2temp,1);

figure(1)
plot(KTAS, CL, '-o')
grid minor
xlabel('Airspeed (kts)');
ylabel('C_L')
hold on
plot(KTAS, CL1, '-*')
plot(KTAS, CL2, '-^')
hold off
legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2')


figure(2)
plot(KTAS, CD, '-o')
grid minor
xlabel('Airspeed (kts)');
ylabel('C_D')
hold on
plot(KTAS, CD1, '-*')
plot(KTAS, CD2, '-^')
hold off
legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2')


figure(3)
plot(KTAS, CT, '-o')
grid minor
xlabel('Airspeed (kts)');
ylabel('C_T')
hold on
plot(KTAS, CT1, '-*')
plot(KTAS, CT2, '-^')
hold off
legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
    'Location','SouthEast')










