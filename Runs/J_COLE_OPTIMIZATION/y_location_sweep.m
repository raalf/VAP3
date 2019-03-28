clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

% cores = 32;
% parpool(cores,'IdleTimeout',800)
home_dir = pwd;

N_chord = 11;
N_dihe = 0;
N_prop_max = 1;
Vars_prop = 3;

max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

stations = 30;
yloc = [linspace(150, 400, stations/2) linspace(410,482, stations/2)];
design_vec = repmat([linspace(75.6, 53, 11) nan nan nan 0 nan],stations,1);
design_vec(:,14) = yloc';

small_prop = 100;
small_prop_rpm = 4153;

big_prop = 160;
big_prop_rpm = 2250;

% % Small prop, inboard down
% z = design_vec;
% z(:,12) = small_prop;
% z(:,13) = small_prop_rpm;
% z(:,16) = 0;
% parfor i = 1:size(z,1)
%     [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
%     
%     sd_out(i).out = out;
%     sd_out(i).ITER = ITER;
%     sd_out(i).ITEROUTP = ITEROUTP;
% end
% 
% % Small prop, inboard up
% z = design_vec;
% z(:,12) = small_prop;
% z(:,13) = small_prop_rpm;
% z(:,16) = 1;
% parfor i = 1:size(z,1)
%     [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
%     
%     su_out(i).out = out;
%     su_out(i).ITER = ITER;
%     su_out(i).ITEROUTP = ITEROUTP;
% end
% 
% % Big prop, inboard down
% z = design_vec;
% z(:,12) = big_prop;
% z(:,13) = big_prop_rpm;
% z(:,16) = 0;
% parfor i = 1:size(z,1)
%     [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
%     
%     bd_out(i).out = out;
%     bd_out(i).ITER = ITER;
%     bd_out(i).ITEROUTP = ITEROUTP;
% end
% 
% % Big prop, inboard up
% z = design_vec;
% z(:,12) = big_prop;
% z(:,13) = big_prop_rpm;
% z(:,16) = 1;
% parfor i = 1:size(z,1)
%     [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
%     
%     bu_out(i).out = out;
%     bu_out(i).ITER = ITER;
%     bu_out(i).ITEROUTP = ITEROUTP;
% end
% 
% save('y_sweep_prop_in_loop.mat')

% load('y_sweep_prop_in_loop.mat')
load('y_sweep_2_20.mat')

hFig24 = figure(24);
clf(24);


tmp = [bu_out.out]'; tmp(tmp > 1e+6) = nan;
plot(yloc./482, smooth(tmp), '--ok');
hold on
tmp = [bd_out.out]'; tmp(tmp > 1e+6) = nan;
plot(yloc./482, smooth(tmp), '-sb');
tmp = [su_out.out]'; tmp(tmp > 1e+6) = nan;
plot(yloc./482, smooth(tmp), '--^m');
tmp = [sd_out.out]'; tmp(tmp > 1e+6) = nan;
plot(yloc./482, smooth(tmp), '-.*r');

% plot(yloc, [bu_out.out]', '--ok');
% hold on
% plot(yloc, tmp, '-sb');
% plot(yloc, [su_out.out]', '--^m');
% plot(yloc, [sd_out.out]', '-.*r');

hold off
legend('Large Propeller - Inboard Up','Large Propeller - Inboard Down','Small Propeller - Inboard Up','Small Propeller - Inboard Down','Location','NorthEast');
grid minor
box on
axis tight

xlabel('Propeller Spanwise Hub Location','FontSize',15);
ylabel('Power (W)','FontSize',15);
