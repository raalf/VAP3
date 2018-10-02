clc
clear

try
    parpool
end

delete optihistory.txt
delete dvhistory.txt

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints_validation;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation_validation(nvars,200,ub,lb,A,b,N_chord,Vars_prop);

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed},nvars,A,b,[],[],lb,ub,[],IntCon,options);
