% clc
clear
warning off

filename = 'inputs/Goland_Wing.vap';

VAP_IN = [];
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