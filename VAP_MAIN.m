clc
clear
warning off

filename = 'inputs/J_COLE_BASELINE_SYM.vap';
% filename = 'inputs/J_COLE_X57_CRUISE_PROP.vap'
% filename = 'inputs/J_COLE_X57_CRUISE_PROP_20_SECTIONS.vap'
 
% seqALPHA = [2:1:15];
% vecCOLLECTIVE = [-50:2:10];
vecCOLLECTIVE = 0;

OUTP = fcnVAP_MAIN(filename, 5, 0);

% for i = 1:length(vecCOLLECTIVE)
%     OUTP(i) = fcnVAP_MAIN(filename, 5, vecCOLLECTIVE(i));
% end

% save('VAP315_STEADY_VISCOUS_FIXED_CRUISE_PROP_20RadialSections')
