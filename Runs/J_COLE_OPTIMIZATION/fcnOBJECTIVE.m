function out = fcnOBJECTIVE(z, N_chord, N_prop_max, Vars_prop)
% clc
% clear
% N_chord = 10;
% N_prop_max = 6;
% Vars_prop = 4;
% z = [100 99 98 97 96 95 94 93 92 91 90, 160 2250, 5 481 0 0, 20 520 0 0, 15 600 0 0, 15 700 0 0, 15 800 0 0, 15 900 0 0];
% z = [75.6000000000000,73.0888888888889,70.5777777777778,68.0666666666667,65.5555555555556,63.0444444444444,60.5333333333333,58.0222222222222,55.5111111111111,53,155,2250,22.6000000000000,482,0,1,32.8215767634855,700,0,1,42.1991701244813,900,0,1,46.8879668049793,1000,0,1,60.9543568464730,1300,0,1,75.0207468879668,1600,0,1];

temp_name = regexprep(tempname('\'), '\', '');
vap_filename = ['aux_files\',temp_name,'.vap'];
copyfile('X57_BLANK.vap', vap_filename);

% Propeller values
rotor.collective = 0;
rotor.axis = [-1 0 0];
rotor.m = 1;
rotor.blades = 3;
rotor.airfoil = 'MH-117';
rotor.diam = z(N_chord + 1)./100;
rotor.rpm = z(N_chord + 2);
airfoil_data = load('airfoils/MH-117.mat');

%% CALCULATED VALUES
for i = 1:N_prop_max
    prop_y(i) = z(N_chord + 2 + (i-1)*Vars_prop + 2);
end

N_prop = sum(prop_y <= 482);

thrust = (13344.66/14)/(2*N_prop); % L/D based on Borer of 14 @ 150 kts, 8000 ft
vinf = 77.2;

%% CREATE PROP IN QMIL
qmil_path = fcnQMILCREATE(temp_name, airfoil_data, rotor.blades, thrust, vinf, rotor.rpm, rotor.diam);
qmil_filename = regexprep(qmil_path, 'aux_files\', '');

% RUN QMIL
qmil_output_filename = ['output_', qmil_filename];
qmil_output_path = ['aux_files\', qmil_output_filename];

exeName = ['aux_files\', qmil_filename, '.exe'];
copyfile('qmil.exe', exeName);

prmpt = sprintf('%s %s %s', exeName, qmil_path, qmil_output_path);
[~,~] = system(prmpt);

delete(exeName);
delete(qmil_path);

%% Build VAP File
geom(:,2) = [0; cumsum(repmat(482/(N_chord-1),N_chord-1, 1))]; % Chord stations
geom(:,1) = (22.6/482).*geom(:,2);
geom(:,5) = (-4/482).*geom(:,2) + 5;
geom(:,3) = geom(:,1).*0;
geom(:,4) = z(1:N_chord)';
geom(:,1:4) = geom(:,1:4)./100; % cm to m

vap3_inputmod_wing(vap_filename, geom)
for i = 1:N_prop
    rotor.hub = z((N_chord + 2 + (i-1)*Vars_prop) + [1:3])./100;
    vap3_inputmod_prop(vap_filename, rotor, qmil_output_path);
end
delete(qmil_output_path);

%% RUN VAP
seqALPHA = 6;
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valSTARTFORCES = 10;
    VAP_IN.valMAXTIME = 10;
    VAP_IN.PREVIEW = 1;
    OUTP(i) = fcnVAP_MAIN(vap_filename, VAP_IN);
end
delete(vap_filename)

%% ANALYZE RESULTS
out = sum(OUTP.vecCP_AVG).*((rotor.rpm/60).^3).*(rotor.diam.^5).*OUTP.valDENSITY;

% end

