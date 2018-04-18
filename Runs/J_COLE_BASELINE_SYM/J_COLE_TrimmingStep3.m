clear
clc
PLOTON = 0;
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

if PLOTON == 1
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
    % saveFig2Latex('JCole_Trim_Rev1_Iter2_CL.pdf', [10, 8])
    
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
    % saveFig2Latex('JCole_Trim_Rev1_Iter2_CD.pdf', [10, 8])
    
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
    % saveFig2Latex('JCole_Trim_Rev1_Iter2_CT.pdf', [10, 8])
end


% delta of targeted CL, CT and calculated CL, CT from 2st iteration
dCL2 = CL1 - CL2;
dCT2 = CT1 - CT2;

% interpolate alpha to maintain steady level flight at VINF 
% using wing+prop data from iter.1 and 2
seqALPHA3 = seqALPHA2*nan;
vecCOLLECTIVE3 = vecCOLLECTIVE2*nan;

if PLOTON == 1
    figure(4)
    clf
    figure(5)
    clf
end

for nn = 1:length(seqALPHA3)
    % New sets of AOA input for 3rd iteration in order to hit the targeted CL
    seqALPHA3(nn) = interp1([CL1(nn) CL2(nn)],...
        [seqALPHA(nn) seqALPHA2(nn)],...
        CL(nn),'linear','extrap');
    
    % New sets of collective pitch input for 3rd iteration in order to hit the targeted CT
    vecCOLLECTIVE3(nn) = interp1([CT1(nn) CT2(nn)],...
        [vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn)],...
        CT(nn),'linear','extrap');
    
    if PLOTON == 1
        figure(4)
        hold on
        plot([seqALPHA(nn) seqALPHA2(nn)],[CL1(nn) CL2(nn)],'-o')
        scatter(seqALPHA3(nn), CL(nn))
        grid minor
        hold off
        
        figure(5)
        hold on
        plot([vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn)],[CT1(nn) CT2(nn)],'-o')
        scatter(vecCOLLECTIVE3(nn), CT(nn))
        grid minor
        hold off
    end
end



% Running
warning off
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE3)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA3(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE3(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    VAP_IN.valSTARTFORCES = 138;
    VAP_IN.valMAXTIME = 160;
    
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);

    fprintf('finished AOA=%.1f\n',seqALPHA(i));
end

save('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps_3rdTrimIter.mat')



