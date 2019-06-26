% clc
clear
warning off

filename = 'inputs/TMotor.vap';

VAP_IN = [];
VAP_IN.RELAX = 0;
VAP_IN.valMAXTIME = 40;
OUTP = fcnVAP_MAIN(filename, VAP_IN);

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