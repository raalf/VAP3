% clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
filename = 'inputs/vap_test.vap'

VAP_IN = [];
OUTP = fcnVAP_MAIN(filename, VAP_IN);

