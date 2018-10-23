% All distance units in cm
% Mixed integer optimization -> inequality constraints only!!
N_chord = 11; % Number of chordwise stations
N_dihe = 0;
N_prop = 1; % Maximum number of propellers
% Constant for all propellers: PROP_DIAM, PROP_RPM
% Unique to each propeller: PROP_XYZ, PROP_DIR
Vars_prop = 3; % [PROP_YZ, PROP_DIR]
Pad_prop = 2 + N_prop*Vars_prop;

Aeq = [];
beq = [];
A = [];
b = [];

%% CHORD LENGTHS
A_area = [];
b_area = [];
A_taper = [];
b_taper = [];

section_length = 482/(N_chord - 1);
area_x57 = ((75.6 + 53)/2)*482; % Half-span area in cm^2

lb_chord = repmat(50, 1, N_chord);
ub_chord = repmat(100, 1, N_chord);
ub_chord(end) = 70;

% Limited negative taper (can increase by 20 cm)
A_taper = -eye(N_chord) + diag(ones(N_chord-1,1),1);
A_taper(end,:) = [];
b_taper = zeros(size(A_taper,1),1) + 20; %cm

% Fixing the area
A_area = ones(2, N_chord);
A_area(:,2:end-1) = A_area(:,2:end-1) + ones(2, N_chord - 2);
A_area(end,:) = A_area(end,:).*-1;

b_area(1,1) = 1.01*(2*area_x57/section_length);
b_area(2,1) = -0.99*(2*area_x57/section_length);

A = [A; padarray([A_taper; A_area], [0 Pad_prop + N_dihe], 0, 'post')];
b = [b; b_taper; b_area];

%% DIHEDRAL
lb_dihe = zeros(1, N_dihe);
ub_dihe = repmat(75, 1, N_dihe);

% % % No anhedral
% A_dihe = eye(N_chord) - diag(ones(N_chord-1,1),1);
% A_dihe(end,:) = [];
% b_dihe = zeros(size(A_dihe,1),1);
% 
% A_dihe =  padarray(A_dihe, [0 Pad_prop], 0, 'post');
% A_dihe =  padarray(A_dihe, [0 N_chord], 0, 'pre');
% A = [A; A_dihe];
% b = [b; b_dihe];

%% PROPELLERS
le_location = 22.6/482;

% xyz of propeller hub, propeller diameter, rpm, rotation direction
lb_prop = [100 2000];
ub_prop = [100 6000];
lb_prop = [lb_prop, repmat([100 -15 0], 1, N_prop)];
ub_prop = [ub_prop, [482 15 1], repmat([2000 15 1], 1, N_prop-1)];

%% COMPILING
lb = [lb_chord, lb_dihe, lb_prop];
ub = [ub_chord, ub_dihe, ub_prop];



