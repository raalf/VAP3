clc
clear

cores = 48;
parpool(cores,'IdleTimeout',800)
home_dir = pwd;

delete opthistory.txt
delete dvhistory.txt

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints_1prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation(nvars,N_prop,191,ub,lb,A,b,N_chord,N_dihe,Vars_prop);
seeds(1,:) = [75.176 71.3419 69.6033 68.1224 64.8988 64.9348 61.9753 59.1098 59.9463 57.9888 52.4117 159.43 2276.93 464.578 -9.70714 0.886147];

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);
