% clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
filename = 'inputs/Goland_Wing.vap';
filename = 'inputs/TMotor.vap';

VAP_IN = [];
VAP_IN.valMAXTIME = 20;
VAP_IN.RELAX = 1;
[OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(filename, VAP_IN);

