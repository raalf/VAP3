clc
clear
warning off

filename = 'inputs/vap_test.vap';
VAP_IN.valMAXTIME = 20
VAP_IN.valSTARTFORCES = 10

OUTP = fcnVAP_MAIN(filename, VAP_IN);


% parfor i = 1:length(seqALPHA)
%     OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
% end
