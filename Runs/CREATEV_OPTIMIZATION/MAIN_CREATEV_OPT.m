clc
clear

cores = 4;
parpool(cores,'IdleTimeout',800)

delete ../../Optimization/opthistory.txt
delete ../../Optimization/dvhistory.txt

warning off

home_dir = pwd;
addpath('../../Flight Dynamics')
addpath('../../')

Constraint_Definition; % Define constraints for all design variables
nvars = length(lb); % Number of design variables
IntCon = [];
[seeds] = fcnPOPCREATION(nvars,200,N_bendstiff,N_torstiff,N_elasticaxis,N_massaxis,ub,lb,A,b); % Generate initial population for optimizer

% Do the thing
options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',seeds,'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_bendstiff, N_torstiff, N_elasticaxis, N_massaxis, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);