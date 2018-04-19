clear
clc

% Import 1st Iteration Data
load('WPforCDo_TS160_22_fixed_iter5.mat')
OUTP5 = OUTP;
clear OUTP;


CL5 = nanmean([OUTP5.vecCL],1);
CT5 = nanmean([OUTP5.vecCT],1);

CD5temp = reshape([OUTP5.vecCD],[],length(OUTP5));
CD5temp(isnan([OUTP5.vecCL])) = nan;
CD5 = nanmean(CD5temp,1);


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
    plot(KTAS, CL5, '-v')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4','Wing + Prop Iter.5')
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
    plot(KTAS, CD5, '-v')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4','Wing + Prop Iter.5')
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
    plot(KTAS, CT5, '-v')
    hold off
    legend('Wing Only','Wing + Prop Iter.1','Wing + Prop Iter.2',...
        'Wing + Prop Iter.3','Wing + Prop Iter.4','Wing + Prop Iter.5'...
        ,'Location','SouthEast')
    %     saveFig2Latex('JCole_Trim_Rev1_Iter3_CT.pdf', [10, 8])
end



% interpolate alpha to maintain steady level flight at VINF
% using wing+prop data from iter.1 and 2
seqALPHA6 = seqALPHA5*nan;
vecCOLLECTIVE6 = vecCOLLECTIVE5*nan;

if PLOTON == 1
    figure(4)
    clf
    figure(5)
    clf
end

for nn = 1:length(seqALPHA5)
    % New sets of AOA input for 3rd iteration in order to hit the targeted CL
    seqALPHA6(nn) = interp1([CL1(nn) CL2(nn) CL3(nn) CL4(nn) CL5(nn)],...
        [seqALPHA(nn) seqALPHA2(nn) seqALPHA3(nn) seqALPHA4(nn) seqALPHA5(nn)],...
        CL(nn),'linear','extrap');
    
    % New sets of collective pitch input for 3rd iteration in order to hit the targeted CT
    vecCOLLECTIVE6(nn) = interp1([CT2(nn) CT3(nn) CT4(nn) CL5(nn)],...
        [vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn) vecCOLLECTIVE4(nn) vecCOLLECTIVE5(nn)],...
        CT(nn),'PCHIP','extrap');
    vecCOLLECTIVE66(nn) = interp1([CT2(nn) CT3(nn) CT4(nn) CT5(nn)],...
        [vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn) vecCOLLECTIVE4(nn) vecCOLLECTIVE5(nn)],...
        CT(nn),'linear','extrap');
    
    if PLOTON == 1
        figure(4)
        hold on
        plot([seqALPHA(nn) seqALPHA2(nn) seqALPHA3(nn) seqALPHA4(nn) seqALPHA5(nn)],...
            [CL1(nn) CL2(nn) CL3(nn) CL4(nn) CL5(nn)],'-o');
        scatter(seqALPHA6(nn), CL(nn),'filled');
        grid minor
        hold off
        
        figure(5)
        hold on
        plot([vecCOLLECTIVE(nn) vecCOLLECTIVE2(nn) vecCOLLECTIVE3(nn) vecCOLLECTIVE4(nn) vecCOLLECTIVE5(nn)],...
            [CT1(nn) CT2(nn) CT3(nn) CT4(nn) CT5(nn)],'-o');
        scatter(vecCOLLECTIVE6(nn), CT(nn),'filled');
        scatter(vecCOLLECTIVE66(nn), CT(nn),'s','filled');
        grid minor
        hold off
    end
end
% 






