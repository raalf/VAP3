clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

% z = [89.7142	80.0991	78.4412	76.4267	73.4695	66.8561	65.5608	63.1486	58.2824	58.1091	50.5951	60.8913	7824.5	98.6364	-1.48375	0.0838214	180.528	-8.13069	0.913337	262.419	-10.4287	0.825817	344.31	1.15027	0.996135	426.202	-12.6547	0.442678];
z = [95.445	92.7811	92.1621	94.3397	87.2277	80.1124	78.1725	66.7732	58.3955	57.8899	55.8115	160	2346.6	435.491	-8.41503	1];
N_chord = 11;
N_dihe = 0;
N_prop_max = 1;
Vars_prop = 3;

max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

[out, ITER, ITEROUTP] = fcnOBJECTIVE(z, N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed);