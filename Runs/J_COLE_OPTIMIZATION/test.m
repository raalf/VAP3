clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

% z = [89.7142	80.0991	78.4412	76.4267	73.4695	66.8561	65.5608	63.1486	58.2824	58.1091	50.5951	60.8913	7824.5	98.6364	-1.48375	0.0838214	180.528	-8.13069	0.913337	262.419	-10.4287	0.825817	344.31	1.15027	0.996135	426.202	-12.6547	0.442678];
% z(1,:) = [95.445	92.7811	92.1621	94.3397	87.2277	80.1124	78.1725	66.7732	58.3955	57.8899	55.8115	160	2346.6	435.491	-8.41503	1];
% z(2,:) = [95.445	92.7811	92.1621	94.3397	87.2277	80.1124	78.1725	66.7732	58.3955	57.8899	55.8115	160	2346.6	482	0	1];
% z(1,:) = [linspace(75.6, 53, 11)	100	4153.52	163.602	-1.00947	0.456371];
% z(2,:) = [linspace(75.6, 53, 11)	100	4153.52	482	-1.00947	1];
% z(3,:) = [linspace(75.6, 53, 11)	100	4153.52	470	-1.00947	1];
z = [98.2576 97.7547 91.4235 91.3494 90.9776 90.6701 86.7798 83.5214 81.0675 69.5841 59.156 62.5367 3985.26 178.649 -11.3943 0.283713 262.186 10.1471 0.901353 345.723 -1.50158 0.988536 ];


N_chord = 11;
N_dihe = 0;
N_prop_max = 3;
Vars_prop = 3;

max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

for i = 1:size(z,1)
        [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
        
        res_out(i).out = out;
        res_out(i).ITER = ITER;
        res_out(i).ITEROUTP = ITEROUTP;

end



% load('matlab.mat')
% 
% for i = 1:2
%         [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
%         
%         res_out(i).out = out;
%         res_out(i).ITER = ITER;
%         res_out(i).ITEROUTP = ITEROUTP;
% 
% end
% 
% hFig278 = figure(278);
% clf(278);
% plot([res_out(1).ITEROUTP(2).OUTP.WING.vecSPANLOC_PROJ; 4.82], [res_out(1).ITEROUTP(2).OUTP.WINGDIST.LDIST; 0], '--ok')
% hold on
% plot([res_out(2).ITEROUTP(2).OUTP.WING.vecSPANLOC_PROJ; 4.82], [res_out(2).ITEROUTP(2).OUTP.WINGDIST.LDIST; 0], '-.rs')
% hold off
% grid minor
% box on
% axis tight
% legend('Unmodified Optimized Case','Propeller at Wingtip','Location','SouthWest')
% xlabel('Spanwise Location (m)','FontSize',15);
% ylabel('Lift (N)','FontSize',15);

% nanmean(res_out(1).ITEROUTP(2).OUTP.vecE)
% nanmean(res_out(2).ITEROUTP(2).OUTP.vecE)