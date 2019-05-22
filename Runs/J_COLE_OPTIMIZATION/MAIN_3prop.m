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

constraints_3prop;
nvars = length(lb);
% IntCon = 1:nvars;
IntCon = [];
[seeds, max_tip_speed, min_tip_speed] = creation(nvars,N_prop,141,ub,lb,A,b,N_chord,N_dihe,Vars_prop);
seeds(1,:) = [76.3609 73.34 71.08 69.358 66.56 64.9754 62.04 59.78 57.5467 55.4369 53 107.184 2980.88 129.922 -9.82695 0.687749 262.291 -1.80909 0.764347 396.786 -5.68468 0.630234]; % 76519;

options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'PopulationSize', 141, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_dihe, N_prop, Vars_prop, max_tip_speed, min_tip_speed, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);
