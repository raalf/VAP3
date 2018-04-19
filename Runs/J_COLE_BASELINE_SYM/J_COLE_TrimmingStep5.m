clear
clc

% Import 1st Iteration Data
load('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps_4thTrimIter.mat')
OUTP4 = OUTP;
clear OUTP;


CL4 = nanmean([OUTP4.vecCL],1);
CT4 = nanmean([OUTP4.vecCT],1);

CD4temp = reshape([OUTP4.vecCD],[],length(OUTP4));
CD4temp(isnan([OUTP4.vecCL])) = nan;
CD4 = nanmean(CD4temp,1);


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
    plot(KTAS, CL4, '-d')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4')
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
    plot(KTAS, CD4, '-d')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4')
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
    plot(KTAS, CT4, '-d')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4','Location','SouthEast')
    %     saveFig2Latex('JCole_Trim_Rev1_Iter3_CT.pdf', [10, 8])
end



% interpolate alpha to maintain steady level flight at VINF
% using wing+prop data from iter.1 and 2
seqALPHA5 = seqALPHA4*nan;
vecCOLLECTIVE5 = vecCOLLECTIVE4*nan;

if PLOTON == 1
    figure(4)
    clf
    figure(5)
    clf
end

for nn = 1:length(seqALPHA5)
    % New sets of AOA input for 3rd iteration in order to hit the targeted CL
    seqALPHA5(nn) = interp1([CL1(nn) CL2(nn) CL3(nn) CL4(nn)],...
        [seqALPHA(nn) seqALPHA2(nn) seqALPHA3(nn) seqALPHA4(nn)],...
        CL(nn),'linear','extrap');
    
    % New sets of collective pitch input for 3rd iteration in order to hit the targeted CT
    vecCOLLECTIVE5(nn) = interp1([CT1(nn) CT2(nn) CT3(nn) CT4(nn)],...
        [vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn) vecCOLLECTIVE4(nn)],...
        CT(nn),'linear','extrap');
    
    if PLOTON == 1
        figure(4)
        hold on
        plot([seqALPHA(nn) seqALPHA2(nn) seqALPHA3(nn) seqALPHA4(nn)],...
            [CL1(nn) CL2(nn) CL3(nn) CL4(nn)],'-o');
        scatter(seqALPHA5(nn), CL(nn),'filled');
        grid minor
        hold off
        
        figure(5)
        hold on
        plot([vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn) vecCOLLECTIVE4(nn)],...
            [CT1(nn) CT2(nn) CT3(nn) CT4(nn)],'-o');
        scatter(vecCOLLECTIVE5(nn), CT(nn),'filled');
        grid minor
        hold off
    end
end


% Running
warning off
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE5)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA5(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE5(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    VAP_IN.valSTARTFORCES = 138;
    VAP_IN.valMAXTIME = 160;
    
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);

    fprintf('finished AOA=%.1f\n',seqALPHA(i));
end

save('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps_5thTrimIter.mat')


















