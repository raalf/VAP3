% function out = fcnOBJECTIVE(z, N_chord, N_prop_max, Vars_prop)
clc
clear
load('temp.mat')
z = [100 99 98 97 96 95 94 93 92 91, 1, -10 480 0 150 2250 1, zeros(1, 6*3)];
% z = [100 99 98 97 96 95 94 93 92 91, 2, -10 480 0 150 2250 1, -10 200 0 110 1400 1, zeros(1, 6*2)];


N_prop = z(N_chord + 1);

temp_name = regexprep(tempname('\'), '\', '');
vap_filename = ['aux_files\',temp_name,'.vap'];
copyfile('X57_WITHOUT_ROTORS.vap', vap_filename);

    %% CALCULATED VALUES
    % Thrust distributed over rotors is proportional to n^2 * d^4
    thrust = (13344.66/14)/(2*N_prop);
    vinf = 80;
    
for i = 1:N_prop
    prop_diam(i) = (z((N_chord + 1 + (i-1)*Vars_prop) + 4))./100;
    prop_rpm(i) = (z((N_chord + 1 + (i-1)*Vars_prop) + 5));
end

    dist = (prop_rpm.^2).*(prop_diam.^4);
    thrust_per_prop = thrust.*(dist./sum(dist));

for i = 1:N_prop
    %% DESIGN VARIABLES (INPUT TO fcnOBJECTIVE)
    rotor.hub = z((N_chord + 1 + (i-1)*Vars_prop) + [1:3])./100;
    rotor.diam = (z((N_chord + 1 + (i-1)*Vars_prop) + 4))./100;
    rotor.rpm = z((N_chord + 1 + (i-1)*Vars_prop) + 5);
    
    rotor.collective = 0;
    rotor.axis = [-1 0 0];
    rotor.m = 1;
    rotor.blades = 3;
    rotor.airfoil = 'MH-117';

    airfoil_data = load('airfoils/MH-117.mat');

    %% CREATE PROP IN QMIL
    qmil_path = fcnQMILCREATE(temp_name, airfoil_data, rotor.blades, thrust_per_prop(i), vinf, rotor.rpm, rotor.diam);
    qmil_filename = regexprep(qmil_path, 'aux_files\', '');

    %% RUN QMIL
    qmil_output_filename = ['output_', qmil_filename];
    qmil_output_path = ['aux_files\', qmil_output_filename];

    exeName = ['aux_files\', qmil_filename, '.exe'];
    copyfile('qmil.exe', exeName);

    prmpt = sprintf('%s %s %s', exeName, qmil_path, qmil_output_path);
    [~,~] = system(prmpt);

    delete(exeName);
    delete(qmil_path);

    %% GENERATE INPUT FILE FOR VAP RUN
    vap3_inputmod(vap_filename, rotor, qmil_output_path);

end

%% RUN VAP
seqALPHA = 12;
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

