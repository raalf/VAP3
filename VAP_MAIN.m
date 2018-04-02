clc
clear
warning off

filename = 'inputs/WinDySIM.vap';

OUTP = fcnVAP_MAIN(filename, 0, 0);

% for i = 1:length(vecCOLLECTIVE)
%     OUTP(i) = fcnVAP_MAIN(filename, 5, vecCOLLECTIVE(i));
% end
