% clc
clear
warning off

filename = 'inputs/VAP_Test.vap';

VAP_IN = [];
VAP_IN.RELAX = 1;
VAP_IN.valMAXTIME = 40;
OUTP = fcnVAP_MAIN(filename, VAP_IN);
