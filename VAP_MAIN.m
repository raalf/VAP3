% clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
filename = 'inputs/Goland_Wing.vap'

VAP_IN = [];
[OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(filename, VAP_IN);

