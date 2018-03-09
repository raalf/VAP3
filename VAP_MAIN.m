clc
clear
warning off

% filename = 'inputs/J_COLE_BASELINE_WING.vap';
filename = 'inputs/J_COLE_X57_CRUISE_PROP.vap'

% seqALPHA = [2:1:15];
vecCOLLECTIVE = [-50:2:10];
% vecCOLLECTIVE = 0;

parfor i = 1:length(vecCOLLECTIVE)
    OUTP(i) = fcnVAP_MAIN(filename, vecCOLLECTIVE(i));
end

save('VAP315ZLL_STEADY_VISCOUS_FIXED_CRUISE_PROP')
