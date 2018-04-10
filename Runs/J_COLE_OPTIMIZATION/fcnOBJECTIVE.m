clc
clear

addpath('../../')
addpath('QMIL')

%% DESIGN VARIABLES (INPUT TO fcnOBJECTIVE)
rotor.rpm = 2250;
rotor.collective = 0;
rotor.dia = 1.524;
rotor.hub = [-0.1 4.82 0];
rotor.axis = [-1 0 0];
rotor.m = 1;
rotor.blades = 3;

%% CREATE PROP IN QMIL


%% GENERATE INPUT FILE FOR VAP RUN
vap_filename = vap3_inputmod(rotor);

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

