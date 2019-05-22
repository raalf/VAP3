clc
clear

cores = 32;
parpool(cores,'IdleTimeout',800)
home_dir = pwd;

if ~exist('aux_files', 'dir')
   mkdir('aux_files');
end

delete opthistory.txt
delete dvhistory.txt

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints_5prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation(nvars,N_prop,160,ub,lb,A,b,N_chord,N_dihe,Vars_prop);

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);
