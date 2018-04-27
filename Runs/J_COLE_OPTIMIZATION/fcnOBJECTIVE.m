function out = fcnOBJECTIVE(z, N_chord, N_prop_max, Vars_prop)

if nargin == 0
    clc
    clear
    N_chord = 10;
    N_prop_max = 6;
    Vars_prop = 4;
    z = [];
end

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

N_prop = sum(prop_y <= 484);

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
    rotor.dir = z((N_chord + 2 + (i-1)*Vars_prop) + 4);
    vap3_inputmod_prop(vap_filename, rotor, qmil_output_path);
end
delete(qmil_output_path);

%% RUN VAP
cd '../../'
seqALPHA = 6;
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valSTARTFORCES = 140;
    VAP_IN.valMAXTIME = 160;
    OUTP(i) = fcnVAP_MAIN(vap_filename, VAP_IN);
    view([90 90]);
end
cd 'Runs/J_COLE_OPTIMIZATION/'
delete(vap_filename)
%% ANALYZE RESULTS
out = sum(2.*OUTP.vecCP_AVG).*((rotor.rpm/60).^3).*(rotor.diam.^5).*OUTP.valDENSITY;

if nargin ~= 0
    fp2 = fopen('opthistory.txt','at');
    fprintf(fp2,'%g ', [out, z]);
    fprintf(fp2,'\n');
    fclose(fp2);
end

end

