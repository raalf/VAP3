clc
clear
warning off

filename = 'inputs/CREATeV.vap';

OUTP = fcnVAP_MAIN(filename, []);


% parfor i = 1:length(seqALPHA)
%     OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
% end
