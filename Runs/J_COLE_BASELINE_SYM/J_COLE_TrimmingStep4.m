clear
clc

% Import 1st Iteration Data
load('VAP32_WPforCDo_TS160_22_fixed_iter3.mat')
OUTP3 = OUTP;
clear OUTP;


CL3 = nanmean([OUTP3.vecCL],1);
CT3 = nanmean([OUTP3.vecCT],1);

CD3temp = reshape([OUTP3.vecCD],[],length(OUTP3));
CD3temp(isnan([OUTP3.vecCL])) = nan;
CD3 = nanmean(CD3temp,1);


PLOTON = 0;

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
    %     saveFig2Latex('JCole_Trim_Rev1_Iter3_CL.pdf', [10, 8])
    
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
    %     saveFig2Latex('JCole_Trim_Rev1_Iter3_CD.pdf', [10, 8])
    
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
    %     saveFig2Latex('JCole_Trim_Rev1_Iter3_CT.pdf', [10, 8])
end



% interpolate alpha to maintain steady level flight at VINF
% using wing+prop data from iter.1 and 2
seqALPHA4 = seqALPHA3*nan;
vecCOLLECTIVE4 = vecCOLLECTIVE3*nan;

if PLOTON == 1
    figure(4)
    clf
    figure(5)
    clf
end

for nn = 1:length(seqALPHA4)
    % New sets of AOA input for 3rd iteration in order to hit the targeted CL
    seqALPHA4(nn) = interp1([CL2(nn) CL3(nn)],...
        [seqALPHA2(nn) seqALPHA3(nn)],...
        CL(nn),'linear','extrap');
    
    % New sets of collective pitch input for 3rd iteration in order to hit the targeted CT
    vecCOLLECTIVE4(nn) = interp1([CT2(nn) CT3(nn)],...
        [vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn)],...
        CT(nn),'linear','extrap');
    
    if PLOTON == 1
        figure(4)
        hold on
        plot([seqALPHA(nn) seqALPHA2(nn) seqALPHA3(nn)],...
            [CL1(nn) CL2(nn) CL3(nn)],'-o');
        scatter(seqALPHA4(nn), CL(nn),'filled');
        grid minor
        hold off
        
        figure(5)
        hold on
        plot([vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn)],...
            [CT1(nn) CT2(nn) CT3(nn)],'-o');
        scatter(vecCOLLECTIVE4(nn), CT(nn),'filled');
        grid minor
        hold off
    end
end



% Running
warning off
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE4)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA4(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE4(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    VAP_IN.valSTARTFORCES = 138;
    VAP_IN.valMAXTIME = 160;
    
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);

    fprintf('finished AOA=%.1f\n',seqALPHA(i));
end

save('VAP32_WPforCDo_TS160_22_fixed_iter4.mat')






