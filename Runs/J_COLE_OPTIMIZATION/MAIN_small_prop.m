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

constraints_small_prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation_validation(nvars,200,ub,lb,A,b,N_chord,Vars_prop);
seeds(1,:) = [75.445 72.566 70.296 67.249 65.547 63.779 60.086 60.699 58.761 60.124 54.47 100 4982.9 463.55 -9.288 0.70461];

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);
