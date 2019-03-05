% clc
clear
warning off

% filename = 'inputs/WIPP.vap';

% VAP_IN = [];
% VAP_IN.valMAXTIME = 1;
% OUTP = fcnVAP_MAIN(filename, VAP_IN);

filename = 'inputs/J_COLE_BASELINE_SYM.vap';

parfor alpha = 1:1:10   
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = alpha
    VAP_IN.valMAXTIME = 160;
    VAP_IN.valSTARTFORCES = VAP_IN.valMAXTIME-20;
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    OUTP(alpha) = fcnVAP_MAIN(filename, VAP_IN);
end

%%
% hFig1 = figure(1);
% clf(1);
% 
% plot([OUTP.vecVEHALPHA]',[OUTP.vecPREQ]','--ok')
% xlabel('Alpha','FontSize',15);
% ylabel('Preq (W)','FontSize',15);
% grid minor
% box on
% axis tight