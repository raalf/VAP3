clc
clear

% All distance units in cm
% Mixed integer optimization -> inequality constraints only!!
N_chord = 11; % Number of chordwise stations
N_prop = 6; % Maximum number of propellers
% Constant for all propellers: PROP_DIAM, PROP_RPM
% Unique to each propeller: PROP_XYZ, PROP_DIR
Vars_prop = 4; % [PROP_XYZ, PROP_DIR]
Pad_prop = 2 + N_prop*Vars_prop;

Aeq = [];
beq = [];
A = [];
b = [];

%% CHORD LENGTHS
A_area = [];
b_area = [];

section_length = 482/(N_chord - 1);
area_x57 = ((75.6 + 53)/2)*482; % Half-span area in cm^2

lb_chord = repmat(50, 1, N_chord);
ub_chord = repmat(100, 1, N_chord);

% No negative taper
A_taper = -eye(N_chord) + diag(ones(N_chord-1,1),1);
A_taper(end,:) = [];
b_taper = zeros(size(A_taper,1),1);

% Fixing the area
A_area = ones(2, N_chord);
A_area(:,2:end-1) = A_area(:,2:end-1) + ones(2, N_chord - 2);
A_area(end,:) = A_area(end,:).*-1;

b_area(1,1) = 1.01*(2*area_x57/section_length);
b_area(2,1) = -0.99*(2*area_x57/section_length);

A = [A; padarray([A_taper; A_area], [0 Pad_prop], 0, 'post')];
b = [b; b_taper; b_area];

%% PROPELLERS
le_location = 22.6/482;

% xyz of propeller hub, propeller diameter, rpm, rotation direction
lb_prop = [75 1300];
ub_prop = [160 2900];
lb_prop = [lb_prop, repmat([-100 100 -30 0], 1, N_prop)];
ub_prop = [ub_prop, [100 482 30 1], repmat([100 2000 30 1], 1, N_prop-1)];

A_prop = [];
b_prop = [];
count = 1;
for i = 1:N_prop
    % Ensuring prop is 10 cm ahead of leading edge at this y-location
    A_prop(count,:) = zeros(1,2+N_prop*Vars_prop);
    A_prop(count, [(i-1)*Vars_prop + 1:(i-1)*Vars_prop + 2]+2) = [1 -le_location];
    b_prop(count, 1) = -20;
    count = count + 1;
    
    % Ensuring prop is within 30 cm of leading edge at this y-location
    A_prop(count,:) = zeros(1,2+N_prop*Vars_prop);
    A_prop(count, [(i-1)*Vars_prop + 1:(i-1)*Vars_prop + 2]+2) = [-1 le_location];
    b_prop(count, 1) = 40;
    count = count + 1;
    
    % Making sure props are spaced along the y-location
    if i < N_prop
        A_prop(count,:) = zeros(1,2+N_prop*Vars_prop);
        A_prop(count, (i-1)*Vars_prop + 4) = 1;
        A_prop(count, 1) = 1;
        A_prop(count, (i)*Vars_prop + 4) = -1;
        b_prop(count, 1) = -20;
        count = count + 1;
    end
end

A = [A; padarray(A_prop, [0 N_chord], 0, 'pre')];
b = [b; b_prop];

%% COMPILING
lb = [lb_chord, lb_prop];
ub = [ub_chord, ub_prop];

%% Seed
Seed = [];
section_location = [0; cumsum(repmat(482/(N_chord-1),N_chord-1, 1))];
x57_chord = section_location.*((53 - 75.6)/482) + 75.6; % Half-span area in cm^2

% X57-ish
y = [482 700 900 1100 1300 1600];
x = le_location.*y - 25;
z = y.*0;
Seed(1,:) = [x57_chord', 155 2250, zeros(1,Vars_prop*N_prop)];
tp = [];
for i = 1:N_prop
   tp = [tp, x(i) y(i) z(i) 1]; 
end
Seed(1,(N_chord + 3):end) = tp;

% Multi-rotor
y = [150 250 350 450 600 700];
x = le_location.*y - 25;
z = y.*0;
direc = [1 0 1 0 1 1];
Seed(2,:) = [x57_chord', 80 3000, zeros(1,Vars_prop*N_prop)];
tp = [];
for i = 1:N_prop
   tp = [tp, x(i) y(i) z(i) direc(i)]; 
end
Seed(2,(N_chord + 3):end) = tp;

% Multi-rotor
y = [200 350 482 900 1300 1600];
x = le_location.*y - 25;
z = y.*0;
direc = [1 0 1 1 1 1];
Seed(3,:) = [x57_chord', 110 2100, zeros(1,Vars_prop*N_prop)];
tp = [];
for i = 1:N_prop
   tp = [tp, x(i) y(i) z(i) direc(i)]; 
end
Seed(3,(N_chord + 3):end) = tp;

% Multi-rotor
y = [250 482 900 1300 1600 1800];
x = le_location.*y - 25;
z = y.*0;
direc = [1 0 1 1 1 1];
Seed(4,:) = [x57_chord', 140 2200, zeros(1,Vars_prop*N_prop)];
tp = [];
for i = 1:N_prop
   tp = [tp, x(i) y(i) z(i) direc(i)]; 
end
Seed(4,(N_chord + 3):end) = tp;


