clc
clear
warning off

cd '../../'
filename = 'Runs/CREATeV/CREATeV.vap';
% VAP_IN.valMAXTIME = 20
% VAP_IN.valSTARTFORCES = 18
% seqALPHA = [2:1:14]
seqALPHA = 20
% OUTP = fcnVAP_MAIN(filename, VAP_IN);
for i = 1:length(seqALPHA)
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valMAXTIME = 40;
    VAP_IN.valSTARTFORCES = 1;
    VAP_IN.RELAX = 1;
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
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