clc
clear
warning off

% filename = 'inputs/J_COLE_BASELINE_WING.vap';
filename = 'inputs/J_COLE_X57_CRUISE_PROP.vap'

% seqALPHA = [2:1:15];
vecCOLLECTIVE = [-50:2:15];
% vecCOLLECTIVE = -67

for i = 1:length(vecCOLLECTIVE)
    OUTP(i) = fcnVAP_MAIN(filename, vecCOLLECTIVE(i));
end

