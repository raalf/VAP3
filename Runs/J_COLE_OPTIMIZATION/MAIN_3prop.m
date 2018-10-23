clc
clear

parpool(32,'IdleTimeout',800)

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints_3prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation(nvars,N_prop,141,ub,lb,A,b,N_chord,N_dihe,Vars_prop);

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'PopulationSize', 141, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed},nvars,A,b,[],[],lb,ub,[],IntCon,options);
