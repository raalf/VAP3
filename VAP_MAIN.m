clc
clear
warning off

filename = 'inputs/J_COLE_BASELINE_WING.vap';

% OUTP = fcnVAP_MAIN(filename, 10, 0);
seqALPHA = 0:2:16;


parfor i = 1:length(seqALPHA)
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
end
