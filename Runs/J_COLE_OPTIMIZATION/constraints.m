clc
clear

% All distance units in cm
% Mixed integer optimization -> inequality constraints only!!
N_chord = 10; % Number of chordwise stations
N_prop = 4; % Maximum number of propellers
Vars_prop = 6; % [PROP_XYZ, PROP_DIAM, PROP_RPM, PROP_DIR]
Pad_prop = 1 + N_prop*Vars_prop;

Aeq = [];
beq = [];
A = [];
b = [];

%% CHORD LENGTHS
A_area = [];
b_area = [];

section_length = 482/N_chord;
area_x57 = ((75.6 + 53)/2)*482; % Half-span area in cm^2

lb_chord = repmat(10, 1, N_chord);
ub_chord = repmat(100, 1, N_chord);

% No negative taper
A_taper = -eye(N_chord) + diag(ones(N_chord-1,1),1);
A_taper(end,:) = [];
b_taper = zeros(size(A_taper,1),1);
% 
% % Fixing the area
% A_area = 0.5*ones(2, N_chord);
% b_area = repmat(area_x57/section_length,2,1);
% A_area(end,:) = A_area(end,:).*-1;
% b_area(end,1) = b_area(end).*-1;

A = [A; padarray([A_taper; A_area], [0 Pad_prop], 0, 'post')];
b = [b; b_taper; b_area];

%% PROPELLERS
le_location = 22.6/482;

% xyz of propeller hub, propeller diameter, rpm, rotation direction
lb_prop = repmat([-30 100 -30 100 1400 0], 1, N_prop);
ub_prop = repmat([20 482 30 200 4000 1], 1, N_prop);

A_prop = [];
b_prop = [];
count = 1;
for i = 1:N_prop
    A_prop(count,:) = zeros(1,N_prop*Vars_prop);
    % Ensuring prop is 10 cm ahead of leading edge at this y-location
    A_prop(count, (i-1)*Vars_prop + 1:(i-1)*Vars_prop + 2) = [1 le_location];
    b_prop(count,1) = -10;
    count = count + 1;
    
    % Prop y-location
    if i < N_prop
        A_prop(count,:) = zeros(1,N_prop*Vars_prop);
        A_prop(count, (i-1)*Vars_prop + 2) = 1;
        A_prop(count, (i-1)*Vars_prop + 4) = 0.5;
        A_prop(count, (i)*Vars_prop + 2) = -1;
        A_prop(count, (i)*Vars_prop + 4) = 0.5;
        b_prop(count,1) = -10;
        count = count + 1;
    end
end

A = [A; padarray(A_prop, [0 N_chord + 1], 0, 'pre')];
b = [b; b_prop];

%% COMPILING
lb = [lb_chord, 4, lb_prop];
ub = [ub_chord, 4, ub_prop];





