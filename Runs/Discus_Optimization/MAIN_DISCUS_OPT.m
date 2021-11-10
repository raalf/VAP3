clc
clear

% cores = 32;
% parpool(cores,'IdleTimeout',800)

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
options = optimoptions('ga', 'Display', 'iter','InitialPopulationMatrix',seeds,'PopulationSize',200,'UseParallel',false,'MaxGenerations',1000, 'StallGenLimit',50,'MutationFcn','mutationadaptfeasible',...
    'CreationFcn','gacreationlinearfeasible','CrossoverFcn',@crossoverintermediate);
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_bendstiff, N_torstiff, N_elasticaxis, N_massaxis, home_dir},nvars,A,b,[],[],lb,ub,[],IntCon,options);

% options = optimoptions('surrogateopt', 'Display', 'iter','InitialPoints',seeds,'MinSurrogatePoints',200,'UseParallel',true,'MaxFunctionEvaluations',6000);
% [x,fval,exitflag,output,trials] = surrogateopt(@fcnOBJECTIVE_SURR,lb,ub,IntCon,A,b,[],[],options);

% options = optimoptions('patternsearch','Display', 'iter', 'UseParallel', true, 'MaxIterations', 6000,'UseCompleteSearch',true,'UseCompletePoll',true,'PollMethod','GSSPositiveBasis2N');
% [x,fval,exitflag,output] = patternsearch({@fcnOBJECTIVE, N_bendstiff, N_torstiff, N_elasticaxis, N_massaxis, home_dir},seeds(1,:),A,b,[],[],lb,ub,IntCon,options);