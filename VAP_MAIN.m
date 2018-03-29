clc
clear
warning off

filename = 'inputs/TMotor-Propeller.vap';

OUTP = fcnVAP_MAIN(filename, 0, 0);
% seqALPHA = 0:1:13;


% parfor i = 1:length(seqALPHA)
%     OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
% end
