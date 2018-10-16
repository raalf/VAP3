clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

N_chord = 11;
N_dihe = 0;
N_prop_max = 1;
Vars_prop = 3;

max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

stations = 15;
yloc = linspace(150,482, stations);
design_vec = repmat([linspace(75.6, 53, 11) nan nan nan 0 nan],stations,1);
design_vec(:,14) = yloc';

small_prop = 100;
small_prop_rpm = 4153;

big_prop = 160;
big_prop_rpm = 2250;

% Small prop, inboard down
z = design_vec;
z(:,12) = small_prop;
z(:,13) = small_prop_rpm;
z(:,16) = 0;
parfor i = 1:size(z,1)
    [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
    
    sd_out(i).out = out;
    sd_out(i).ITER = ITER;
    sd_out(i).ITEROUTP = ITEROUTP;
end

% Small prop, inboard up
z = design_vec;
z(:,12) = small_prop;
z(:,13) = small_prop_rpm;
z(:,16) = 1;
parfor i = 1:size(z,1)
    [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
    
    su_out(i).out = out;
    su_out(i).ITER = ITER;
    su_out(i).ITEROUTP = ITEROUTP;
end

% Big prop, inboard down
z = design_vec;
z(:,12) = big_prop;
z(:,13) = big_prop_rpm;
z(:,16) = 0;
parfor i = 1:size(z,1)
    [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
    
    bd_out(i).out = out;
    bd_out(i).ITER = ITER;
    bd_out(i).ITEROUTP = ITEROUTP;
end

% Big prop, inboard up
z = design_vec;
z(:,12) = big_prop;
z(:,13) = big_prop_rpm;
z(:,16) = 1;
parfor i = 1:size(z,1)
    [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);
    
    bu_out(i).out = out;
    bu_out(i).ITER = ITER;
    bu_out(i).ITEROUTP = ITEROUTP;
end

save('y_sweep.mat')