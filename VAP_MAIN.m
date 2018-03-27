clc
clear
warning off

filename = 'inputs/J_COLE_BASELINE_WING.vap';

OUTP = fcnVAP_MAIN(filename, 5, 0);

% for i = 1:length(vecCOLLECTIVE)
%     OUTP(i) = fcnVAP_MAIN(filename, 5, vecCOLLECTIVE(i));
% end
