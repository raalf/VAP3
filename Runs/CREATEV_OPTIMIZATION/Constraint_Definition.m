N_bendstiff = 15;
N_torstiff = 15;
N_elasticaxis = 15;
N_massaxis = 15;

A = [];
b = [];
lb = [];
ub = [];

%% Bending stiffness constraints
% Bending stiffness from one station to the next cannot change by more than
% 5% of the previous value
tempEI1 = diag(-1.05*ones(N_bendstiff,1),0);
tempEI2 = diag(ones(N_bendstiff-1,1),1);

A_bendstiff = tempEI1 + tempEI2;
b_bendstiff = zeros(N_bendstiff,1);

% Bending stiffness must be greater than 5000 Nm2
lb_bendstiff = 5000*ones(N_bendstiff,1);
ub_bendstiff = Inf*ones(N_bendstiff,1);

A = [A; A_bendstiff];
b = [b; b_bendstiff];

lb = [lb; lb_bendstiff];
ub = [ub; ub_bendstiff];

%% Torsional stiffness constraints
% Bending stiffness from one station to the next cannot change by more than
% 5% of the previous value
tempGJ1 = diag(-1.05*ones(N_torstiff,1),0);
tempGJ2 = diag(ones(N_torstiff-1,1),1);

A_torstiff = tempGJ1 + tempGJ2;
b_torstiff = zeros(N_torstiff,1);

% Torsional stiffness must be greater than 5000 Nm2
lb_torstiff = 5000*ones(N_torstiff,1);
ub_torstiff = Inf*ones(N_torstiff,1);

A = [padarray(A, [0 N_torstiff], 0, 'post'); padarray(A_torstiff,[0 N_bendstiff], 0, 'pre')];
b = [b; b_torstiff];

lb = [lb; lb_torstiff];
ub = [ub; ub_torstiff];

%% Elastic axis constraints
% Chordwise location from one station to the next cannot change by more
% than 2.5% of the previous value
tempEA1 = diag(-1.025*ones(N_elasticaxis,1),0);
tempEA2 = diag(ones(N_elasticaxis-1,1),1);

A_elasticaxis = tempEA1 + tempEA2;
b_elasticaxis = zeros(N_elasticaxis,1);

% Elastic axis must be between 0.25c and 0.75c
lb_elasticaxis = 0.25*ones(N_elasticaxis,1);
ub_elasticaxis = 0.75*ones(N_elasticaxis,1);

A = [padarray(A, [0 N_elasticaxis], 0, 'post'); padarray(A_elasticaxis,[0 N_torstiff+N_bendstiff], 0, 'pre')];
b = [b; b_elasticaxis];

lb = [lb; lb_elasticaxis];
ub = [ub; ub_elasticaxis];

%% Mass axis constraints
% Chordwise location from one station to the next cannot change by more
% than 2.5% of the previous value
tempCG1 = diag(-1.025*ones(N_massaxis,1),0);
tempCG2 = diag(ones(N_massaxis-1,1),1);

A_massaxis = tempCG1 + tempCG2;
b_massaxis = zeros(N_massaxis,1);

% Mass axis must be between 0.25c and 0.75c
lb_massaxis = 0.25*ones(N_massaxis,1);
ub_massaxis = 0.75*ones(N_massaxis,1);

A = [padarray(A, [0 N_massaxis], 0, 'post'); padarray(A_massaxis,[0 N_elasticaxis+N_torstiff+N_bendstiff], 0, 'pre')];
b = [b; b_massaxis];

lb = [lb; lb_massaxis];
ub = [ub; ub_massaxis];

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