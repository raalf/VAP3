clc
clear

delete optihistory.txt

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

warning off

constraints;
nvars = length(lb);
IntCon = 1:nvars;

TolCon_Data = 1e-6; 
TolFun_Data = 1e-08;

load('seeds.mat');

options = optimoptions('ga','InitialPopulationMatrix', Seeds, 'Display', 'iter', 'UseVectorized', false, 'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'PopulationSize', 120);
[x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, N_chord, N_prop, Vars_prop},nvars,A,b,[],[],lb,ub,[],IntCon,options);