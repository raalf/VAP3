clc
clear

seqALPHA = [-4:1:12];

seqALPHA = 13

filename = 'inputs/CREATeV.vap';
for i = 1:length(seqALPHA)
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end
% 
% save('matlab2.mat')
% load('matlab2.mat')
% 
% hFig20 = figure(20);
% clf(20);
% plot([OUTP.vecVEHALPHA],[OUTP.vecCLv_AVG],'-ok')
% xlabel('Vehicle Angle of Attack','FontSize',15);
% ylabel('Lift Coefficient','FontSize',15);
% grid minor
% box on
% 
% hFig21 = figure(21);
% clf(21);
% plot([OUTP.vecCD_AVG], [OUTP.vecCLv_AVG],'-ok')
% xlabel('Drag Coefficient','FontSize',15);
% ylabel('Lift Coefficient','FontSize',15);
% grid minor
% box on

load('airfoils/RAALF-B8.mat')

linestyles = {'--';'-.';'-';':'};
markers = {'o';'x';'s';'^';'*';'d';'v';'>';'<';'p';'h'};
colors = {'k';'b';'r';'m';'c';'g'};

hFig23 = figure(23);
clf(23);
hold on
for i = 1:size(pol,3)
   
    plot(pol(:,1,i), pol(:,2,i), [linestyles{1+mod(i,4)}, markers{1+mod(i,11)}, colors{1+mod(i,6)}])
    
end
hold off
grid minor
box on
