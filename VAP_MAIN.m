clc
clear
warning off

filename = 'inputs/J_COLE_BASELINE_WING.vap';

OUTP = fcnVAP_MAIN(filename, 5, 0);
% seqALPHA = 0:1:13;


% parfor i = 1:length(seqALPHA)
%     OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
% end
