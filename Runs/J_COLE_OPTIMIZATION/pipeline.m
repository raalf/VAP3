function out = fcnOBJECTIVE(z)

%% DESIGN VARIABLES (INPUT TO fcnOBJECTIVE)
rotor.rpm = 2250;
rotor.collective = 0;
rotor.hub = [-0.1 4.82 0];
rotor.axis = [-1 0 0];
rotor.m = 1;
rotor.rpm = 2000;
rotor.diam = 1.524;
rotor.blades = 3;
rotor.airfoil = 'MH-117';

airfoil_data = load('airfoils/MH-117.mat')
%% CALCULATED VALUES
thrust = 600
vinf = 80

%% CREATE PROP IN QMIL
qmil_path = fcnQMILCREATE(airfoil_data, rotor.blades, thrust, vinf, rotor.rpm, rotor.diam);
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
vap_filename = vap3_inputmod(rotor, qmil_output_path);

%% RUN VAP
seqALPHA = 12;
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valSTARTFORCES = 29;
    VAP_IN.valMAXTIME = 30;
        
    OUTP(i) = fcnVAP_MAIN(vap_filename, VAP_IN);
end
delete(vap_filename)

%% ANALYZE RESULTS

