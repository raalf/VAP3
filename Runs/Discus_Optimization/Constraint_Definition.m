N_bendstiff = 19;
N_torstiff = 19;
N_elasticaxis = 19;
N_massaxis = 19;
N_lmass = 19;

A = [];
b = [];
lb = [];
ub = [];

%% Bending stiffness constraints
% Bending stiffness from one station to the next cannot change by more than
% 5% of the previous value

% Setting up upper bound of relative change constraints
tempEI1 = diag(-1.05*ones(N_bendstiff,1),0);
tempEI2 = diag(ones(N_bendstiff-1,1),1);

A_bendstiff = tempEI1 + tempEI2;
b_bendstiff = zeros(N_bendstiff,1);

A_bendstiff(end,end) = 0;

A = [A; A_bendstiff];
b = [b; b_bendstiff];

% Setting up lower bound of relative change constraints
tempEI1 = diag(0.95*ones(N_bendstiff,1),0);
tempEI2 = diag(-1*ones(N_bendstiff-1,1),1);

A_bendstiff1 = tempEI1 + tempEI2;
b_bendstiff1 = zeros(N_bendstiff,1);

A_bendstiff1(end,end) = 0;

% Bending stiffness must be greater than 20000 Nm2
lb_bendstiff = 20000*ones(N_bendstiff,1);
ub_bendstiff = 1e6*ones(N_bendstiff,1);

lb = [lb; lb_bendstiff];
ub = [ub; ub_bendstiff];

A = [A; A_bendstiff1];
b = [b; b_bendstiff1];

%% Torsional stiffness constraints
% Bending stiffness from one station to the next cannot change by more than
% 5% of the previous value

% Setting up upper bound of relative change constraints
tempGJ1 = diag(-1.05*ones(N_torstiff,1),0);
tempGJ2 = diag(ones(N_torstiff-1,1),1);

A_torstiff = tempGJ1 + tempGJ2;
b_torstiff = zeros(N_torstiff,1);

A_torstiff(end,end) = 0;

A = [padarray(A, [0 N_torstiff], 0, 'post'); padarray(A_torstiff,[0 N_bendstiff], 0, 'pre')];
b = [b; b_torstiff];

% Setting up lower bound of relative change constraints
tempGJ1 = diag(0.95*ones(N_torstiff,1),0);
tempGJ2 = diag(-1*ones(N_torstiff-1,1),1);

A_torstiff1 = tempGJ1 + tempGJ2;
b_torstiff1 = zeros(N_torstiff,1);

A_torstiff1(end,end) = 0;

% Torsional stiffness must be greater than 20000 Nm2
lb_torstiff = 20000*ones(N_torstiff,1);
ub_torstiff = 1e6*ones(N_torstiff,1);

lb = [lb; lb_torstiff];
ub = [ub; ub_torstiff];

A = [A; padarray(A_torstiff1, [0 N_bendstiff], 0, 'pre')];
b = [b; b_torstiff1];

%% Elastic axis constraints
% Chordwise location from one station to the next cannot change by more
% than 2.5% of the previous value

% Setting up upper bound of relative change constraints
tempEA1 = diag(-1.025*ones(N_elasticaxis,1),0);
tempEA2 = diag(ones(N_elasticaxis-1,1),1);

A_elasticaxis = tempEA1 + tempEA2;
b_elasticaxis = zeros(N_elasticaxis,1);

A_elasticaxis(end,end) = 0;

A = [padarray(A, [0 N_elasticaxis], 0, 'post'); padarray(A_elasticaxis,[0 N_torstiff+N_bendstiff], 0, 'pre')];
b = [b; b_elasticaxis];

% Setting up lower bound of relative change constraints
tempEA1 = diag(0.975*ones(N_elasticaxis,1),0);
tempEA2 = diag(-1*ones(N_elasticaxis-1,1),1);

A_elasticaxis1 = tempEA1 + tempEA2;
b_elasticaxis1 = zeros(N_elasticaxis,1);

A_elasticaxis1(end,end) = 0;

% Elastic axis must be between 0.20c and 0.75c
lb_elasticaxis = 0.20*ones(N_elasticaxis,1);
ub_elasticaxis = 0.75*ones(N_elasticaxis,1);

lb = [lb; lb_elasticaxis];
ub = [ub; ub_elasticaxis];

A = [A; padarray(A_elasticaxis1, [0 N_bendstiff+N_torstiff], 0, 'pre')];
b = [b; b_elasticaxis1];

%% Mass axis constraints
% Chordwise location from one station to the next cannot change by more
% than 2.5% of the previous value

% Setting up upper bound of relative change constraints
tempCG1 = diag(-1.025*ones(N_massaxis,1),0);
tempCG2 = diag(ones(N_massaxis-1,1),1);

A_massaxis = tempCG1 + tempCG2;
b_massaxis = zeros(N_massaxis,1);

A_massaxis(end,end) = 0;

A = [padarray(A, [0 N_massaxis], 0, 'post'); padarray(A_massaxis,[0 N_elasticaxis+N_torstiff+N_bendstiff], 0, 'pre')];
b = [b; b_massaxis];

% Setting up lower bound of relative change constraints
tempCG1 = diag(0.975*ones(N_massaxis,1),0);
tempCG2 = diag(-1*ones(N_massaxis-1,1),1);

A_massaxis1 = tempCG1 + tempCG2;
b_massaxis1 = zeros(N_massaxis,1);

A_massaxis1(end,end) = 0;

% Mass axis must be between 0.20c and 0.75c
lb_massaxis = 0.20*ones(N_massaxis,1);
ub_massaxis = 0.75*ones(N_massaxis,1);

lb = [lb; lb_massaxis];
ub = [ub; ub_massaxis];

A = [A; padarray(A_massaxis1, [0 N_bendstiff+N_torstiff+N_elasticaxis], 0, 'pre')];
b = [b; b_massaxis1];

%% Constraints between mass and elastic axis
% Constrain mass axis to be in front of elastic axis
A_masselastic = diag(-1*ones(N_elasticaxis,1),0);
A_masselastic = [A_masselastic,diag(ones(N_massaxis,1),0)];

b_masselastic = zeros(N_elasticaxis,1);

A = [A; padarray(A_masselastic, [0 N_bendstiff+N_torstiff], 0, 'pre')];
b = [b; b_masselastic];

% Constrain mass axis and elastic axis to be within 0.3c of each other
A_masselastic1 = diag(ones(N_elasticaxis,1),0);
A_masselastic1 = [A_masselastic1,diag(-1*ones(N_massaxis,1),0)];

b_masselastic1 = 0.3*ones(N_elasticaxis,1);

A = [A; padarray(A_masselastic1, [0 N_bendstiff+N_torstiff], 0, 'pre')];
b = [b; b_masselastic1];

%% Mass constraints
% Chordwise location from one station to the next cannot change by more
% than 2.5% of the previous value

% Setting up upper bound of relative change constraints
% tempLM1 = diag(-1.025*ones(N_lmass,1),0);
% tempLM2 = diag(ones(N_lmass-1,1),1);
% 
% A_lmass = tempLM1 + tempLM2;
% b_lmass = zeros(N_lmass,1);
% 
% A_lmass(end,end) = 0;
% 
% A = [padarray(A, [0 N_lmass], 0, 'post'); padarray(A_lmass,[0 N_elasticaxis+N_torstiff+N_bendstiff+N_massaxis], 0, 'pre')];
% b = [b; b_massaxis];
% 
% % Setting up lower bound of relative change constraints
% tempLM1 = diag(0.975*ones(N_massaxis,1),0);
% tempLM2 = diag(-1*ones(N_massaxis-1,1),1);
% 
% A_lmass1 = tempLM1 + tempLM2;
% b_lmass1 = zeros(N_lmass,1);
% 
% A_lmass1(end,end) = 0;
% 
% % Mass axis must be between 0.20c and 0.75c
% lb_lmass = 0.30*ones(N_lmass,1);
% ub_lmass = 2.5*ones(N_lmass,1);
% 
% lb = [lb; lb_lmass];
% ub = [ub; ub_lmass];
% 
% A = [A; padarray(A_lmass1, [0 N_bendstiff+N_torstiff+N_elasticaxis+N_massaxis], 0, 'pre')];
% b = [b; b_lmass1];