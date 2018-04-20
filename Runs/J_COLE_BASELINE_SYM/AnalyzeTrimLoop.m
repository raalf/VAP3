% load('ITER_0.mat');
load('JCOLE_TRIM_LOOP.mat');
% ITER.maxIter = 5;
% ITER.numCase = 5;

figure(1)
clf
plot(KTAS, CL, '-o')
hold on
plot(repmat(KTAS,ITER.maxIter,1)', ITER.CL', ':x')
hold off
grid minor
legend('Borer Data','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3')
xlabel('Airspeed (kts)');
ylabel('C_L')

figure(2)
clf
plot(KTAS, CT, '-o')
hold on
plot(repmat(KTAS,ITER.maxIter,1)', ITER.CT', ':x')
hold off
grid minor
legend('Borer Data','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3', 'Location', 'northwest')
xlabel('Airspeed (kts)');
ylabel('C_T')

figure(3)
clf
plot(KTAS, CD, '-o')
hold on
plot(repmat(KTAS,ITER.maxIter,1)', ITER.CD', ':x')
hold off
grid minor
xlabel('Airspeed (kts)');
ylabel('C_D')
legend('Borer Data','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3')
ylim([0 0.08]);

% figure(4)
% plot(ITER.AOA, ITER.CL)
% hold on
% scatter(ITER.AOA(:), ITER.CL(:), 20, ITER.Iteration(:), 'filled')
% % scatter(seqALPHA, CL, 'd')
% hold off
% grid minor
% xlabel('Alpha, deg')
% ylabel('C_L')
% 
% figure(5)
% plot(ITER.CLTV, ITER.CT)
% hold on
% scatter(ITER.CLTV(:), ITER.CT(:), 20, ITER.Iteration(:), 'filled')
% % scatter(vecCOLLECTIVE, CT, 'd')
% hold off
% grid minor
% xlabel('Collective Pitch, deg')
% ylabel('C_T')
% 
% 
figure(6)
plot(KTAS, CD-ITER.CD(3,:), '-om')
hold on
plot(repmat(KTAS,ITER.maxIter,1)', (repmat(CD,ITER.maxIter,1)-ITER.CD)', ':x')
hold off
ylim([0 0.05]);
grid minor
xlabel('Airspeed (kts)');
ylabel('C_D_o')
% 
% 
% figure(7)
% plot(ITER.CD', ITER.CL', '-')
% grid minor
% % xlabel('Airspeed (kts)');
% % ylabel('C_D_o')


CL_Delta = (ITER.CL-repmat(CL,ITER.maxIter,1))./repmat(CL,ITER.maxIter,1)
CT_Delta = (ITER.CT-repmat(CT,ITER.maxIter,1))./repmat(CT,ITER.maxIter,1)

CDo = CD-ITER.CD(3,:)

