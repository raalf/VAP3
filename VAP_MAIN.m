clc
clear
warning off

filename = 'inputs/WinDySIM.vap';

OUTP = fcnVAP_MAIN(filename, 0, 0);


% parfor i = 1:length(seqALPHA)
%     OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
% end
