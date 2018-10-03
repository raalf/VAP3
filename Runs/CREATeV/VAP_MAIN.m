clc
clear
warning off

cd G:\GIT\VAP3\Runs\CREATeV

%%
RELAX = 0
valMAXTIME = 40;
valSTARTFORCES = 40;
seqALPHA = [6:1:15];

%%
% a = dir;
% for i = 1:size(a,1)
%     if contains(a(i).name,'Design_')
%
%         design_num = cellfun(@str2num,regexp(a(i).name, '[0-9]','match'));
%
%         cd '../../'
%         filename = ['Runs/CREATeV/', a(i).name];
%         parfor j = 1:length(seqALPHA)
%             VAP_IN = [];
%             VAP_IN.vecVEHALPHA = seqALPHA(j);
%             VAP_IN.valMAXTIME = valMAXTIME;
%             VAP_IN.valSTARTFORCES = valSTARTFORCES;
%             VAP_IN.RELAX = RELAX;
%             OUTP(j) = fcnVAP_MAIN(filename, VAP_IN);
%         end
%
%         Design(design_num).OUTP = OUTP;
%
%         cd 'Runs/CREATeV';
%
%     end
% end
% 
cd '../../'
filename = ['Runs/CREATeV/Design_', num2str(6),'.vap'];
parfor j = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(j);
    VAP_IN.valMAXTIME = valMAXTIME;
    VAP_IN.valSTARTFORCES = valSTARTFORCES;
    VAP_IN.RELAX = RELAX;
    OUTP(j) = fcnVAP_MAIN(filename, VAP_IN);
end
cd 'Runs/CREATeV';

%%
hFig1 = figure(1);
clf(1);

linestyles = {'--ok', '-.rs', '-b^', '-m+'};

hold on
load('design_5_fixed.mat');
plot([OUTP.vecVINF]',[OUTP.vecPREQ]',linestyles{1})
load('design_6_fixed.mat');
plot([OUTP.vecVINF]',[OUTP.vecPREQ]',linestyles{2})
hold off

legend('With Dihedral and Twist','Without Dihedral and Twist','Location','NorthWest')

%
% hold on
% for i = 1:size(Design,2)
%     plot([Design(i).OUTP.vecVINF]',[Design(i).OUTP.vecPREQ]',linestyles{i})
% end
% hold off
xlabel('Velocity (m/s)','FontSize',15);
ylabel('Preq (W)','FontSize',15);
grid minor
box on
axis tight














