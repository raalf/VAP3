clc
clear
warning off

filename = 'inputs/WIPP_PROP.vap';
valAZNUM = 80;
VAP_IN = [];
VAP_IN.vecCOLLECTIVE = 0;
VAP_IN.valMAXTIME = 280;
VAP_IN.valSTARTFORCES = 200;
VAP_IN.vecVEHALPHA = 0;
rps = 4705/60;
VAP_IN.valDELTIME = 1./(valAZNUM.*rps);
OUTP = fcnVAP_MAIN(filename, VAP_IN);
