% clc
% clear
% 
% addpath('../../')
% addpath('../../airfoils')
% addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
% addpath('./../../Runs/J_COLE_OPTIMIZATION/')
% 
% home_dir = pwd;
% 
% N_chord = 11;
% N_dihe = 0;
% N_prop_max = 1;
% Vars_prop = 3;
% 
% max_tip_speed = 277.7821;
% min_tip_speed = 165.3465;
% 
% stations = 10;
% yloc = [linspace(100, 400, stations)];
% yloc = sort([yloc 140.628]);
% stations = size(yloc,2);
% design_vec = repmat([75.6	73.8121	71.5716	69.7713	66.8002	64.6841	62.5955	60.5883	57.52	55.4883	53	100	3622.11	nan	5.64315	0.9],stations,1);
% design_vec(:,14) = yloc';
% 
% % Small prop, inboard down
% z = design_vec;
% for i = 1:size(z,1)
%     [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
%     
%     sd_out(i).out = out;
%     sd_out(i).ITER = ITER;
%     sd_out(i).ITEROUTP = ITEROUTP;
% end
clc
clear

load('prop_sweep_problem.mat')

hFig25 = figure(25);
clf(25);

tmp = [sd_out(1:end).out]';
plot(yloc(1:end)./482, smooth(tmp), '-sk');

load('prop_sweep_other_way.mat')
tmp = [sd_out(2:end).out]';
hold on
plot(yloc(2:end)./482, smooth(tmp), '--b^');
hold off

legend('Inboard Down','Inboard Up','Location','SouthEast')

grid minor
box on
axis tight

xlabel('Propeller Spanwise Hub Location','FontSize',15);
ylabel('Power (W)','FontSize',15);

