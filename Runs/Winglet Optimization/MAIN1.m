clc
clear

clc
clear

% cores = 32;
% parpool(cores,'IdleTimeout',800)
home_dir = pwd;

delete opthistory1.txt
delete dvhistory1.txt

lb = [...
    0.01 0.15 0 ... % Lop location, f trans chord, f trans twist
      0.30 6.8 0.38 0.1 -3 ... % 
      0.30 6.9 0.38 0.05 -5 ...
    ];

ub = [ ...
    0.7 0.37 5 ... 
      0.50 7.49 0.9 0.35 5 ...
      0.50 7.50 0.9 0.35 5 ...
    ];


A = [... 
    -1 0 0 ... % Front tip going outward (inboard)
        0 -1 0 0 0 ...
        0 0 0 0 0; ...
        
    0 0 0 ... % Front tip going outward (outboard)
        0 1 0 0 0 ...
        0 -1 0 0 0; ...            
        
    0 -1 0 ... % Taper front inner
        0 0 0 1 0 ...
        0 0 0 0 0; ...  
        
    0 0 0 ... % Taper front outer
        0 0 0 -1 0 ...
        0 0 0 1 0; ...  
      ];

b = [...
    -7.65; % front tip going outward (inboard)
    -0.03; % front tip going outward (outboard)
    0;
    0;
    ];

Aeq = [];
beq = [];

nvars = 13;
TolCon_Data = 1e-6; 
TolFun_Data = 1e-08;

Seed = [];

options = gaoptimset;
options = gaoptimset(options,'TolFun', TolFun_Data);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'InitialPopulation', Seed);
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv @gaplotexpectation @gaplotscorediversity @gaplotstopping });
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'UseParallel', 1 )
options = gaoptimset(options,'Generations',1000,'StallGenLimit', 50);
[x,fval,exitflag,output,population,score] = gamultiobj({@fcnOBJECTIVE1, home_dir},nvars,A,b,Aeq,beq,lb,ub,[],options);
