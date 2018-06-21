clc
clear

try
    parpool
end


delete optihistory.txt

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints_1prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = []
seeds = [];

options = optimoptions('ga', 'Display', 'iter', 'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'PopulationSize', 120, 'MutationFcn','mutationadaptfeasible', 'CrossoverFraction', 0.5);
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_prop, Vars_prop},nvars,A,b,[],[],lb,ub,[],IntCon,options);
