clc
clear

cores = 32;
parpool(cores,'IdleTimeout',800)

delete ../../Optimization/opthistory_cosine.txt
delete ../../Optimization/dvhistory_cosine.txt

warning off

home_dir = pwd;
addpath('../../Flight Dynamics')
addpath('../../')

Constraint_Definition; % Define constraints for all design variables
nvars = length(lb); % Number of design variables
IntCon = [];
[seeds] = fcnPOPCREATION(nvars,200,N_bendstiff,N_torstiff,N_elasticaxis,N_massaxis,ub,lb,A,b); % Generate initial population for optimizer

% Do the thing
options = optimoptions('ga', 'Display', 'iter','InitialPopulationMatrix',seeds,'PopulationSize',200,'UseParallel',true,'MaxGenerations',1000, 'StallGenLimit',50,'MutationFcn','mutationadaptfeasible',...
    'CreationFcn','gacreationlinearfeasible','CrossoverFcn',@crossoverintermediate);
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE_COSINE, N_bendstiff, N_torstiff, N_elasticaxis, N_massaxis, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);
