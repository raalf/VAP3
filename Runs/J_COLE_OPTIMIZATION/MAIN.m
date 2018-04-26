clc
clear

addpath('../../')
addpath('../../airfoils')

constraints;
nvars = length(lb);
IntCon = 1:nvars;

TolCon_Data = 1e-6; 
TolFun_Data = 1e-08;

Seed = [];

options = gaoptimset;
options = gaoptimset(options,'TolFun', TolFun_Data);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'InitialPopulation', Seed);
options = gaoptimset(options,'PlotFcns', {  @gaplotscorediversity @gaplotstopping @gaplotgenealogy @gaplotscores @gaplotdistance @gaplotselection});
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'populationsize', 200); 
options = gaoptimset(options,'UseParallel', 0);
options = gaoptimset(options,'Generations',1000,'StallGenLimit', 50);
[x,fval,exitflag,output,population,score] = gamultiobj({@fcnOBJECTIVE, N_chord, N_prop, Vars_prop},nvars,A,b,Aeq,beq,lb,ub,[],options);