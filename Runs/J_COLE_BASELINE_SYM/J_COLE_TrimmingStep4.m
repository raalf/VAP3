clear
clc

% Import 1st Iteration Data
load('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps_3rdTrimIter.mat')
OUTP3 = OUTP;
clear OUTP;


CL3 = nanmean([OUTP3.vecCL],1);
CT3 = nanmean([OUTP3.vecCT],1);

CD3temp = reshape([OUTP3.vecCD],[],length(OUTP3));
CD3temp(isnan([OUTP3.vecCL])) = nan;
CD3 = nanmean(CD3temp,1);


PLOTON = 1;
if PLOTON == 1
    figure(1)
    plot(KTAS, CL, '-o')
    grid minor
    xlabel('Airspeed (kts)');
    ylabel('C_L')
    hold on
    plot(KTAS, CL1, '-*')
    plot(KTAS, CL2, '-^')
    plot(KTAS, CL3, '-s')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2','Wing + Prop Iter.3')
    saveFig2Latex('JCole_Trim_Rev1_Iter3_CL.pdf', [10, 8])
    
    figure(2)
    plot(KTAS, CD, '-o')
    grid minor
    xlabel('Airspeed (kts)');
    ylabel('C_D')
    hold on
    plot(KTAS, CD1, '-*')
    plot(KTAS, CD2, '-^')
    plot(KTAS, CD3, '-s')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2','Wing + Prop Iter.3')
    saveFig2Latex('JCole_Trim_Rev1_Iter3_CD.pdf', [10, 8])
    
    figure(3)
    plot(KTAS, CT, '-o')
    grid minor
    xlabel('Airspeed (kts)');
    ylabel('C_T')
    hold on
    plot(KTAS, CT1, '-*')
    plot(KTAS, CT2, '-^')
    plot(KTAS, CT3, '-s')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2','Wing + Prop Iter.3',...
        'Location','SouthEast')
    saveFig2Latex('JCole_Trim_Rev1_Iter3_CT.pdf', [10, 8])
end



