clc
clear

cd G:\VAP3\Runs\J_COLE_OPTIMIZATION
addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

% z = [75.6000000000000,73.8121000000000,71.5716000000000,69.7713000000000,66.8002000000000,64.6841000000000,62.5955000000000,60.5883000000000,57.5200000000000,55.4883000000000,53,100,3622.11000000000,366.666666666667,5.64315000000000,0.282642000000000];
z = [75.6000000000000,73.8121000000000,71.5716000000000,69.7713000000000,66.8002000000000,64.6841000000000,62.5955000000000,60.5883000000000,57.5200000000000,55.4883000000000,53,100,3622.11000000000,100,5.64315000000000,0.900000000000000]

N_chord = 11;
N_dihe = 0;
N_prop_max = 1;
Vars_prop = 3;

% clear
% load('z_test.mat')

max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

home_dir = pwd;
for i = 1:size(z,1)
        [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);

        res_out(i).out = out;
        res_out(i).ITER = ITER;
        res_out(i).ITEROUTP = ITEROUTP;

end




