clc
clear
warning off

filename = 'inputs/J_COLE_BASELINE_WING.vap';
% filename = 'inputs/J_COLE_X57_CRUISE_PROP.vap'

seqALPHA = [2:1:15];

parfor i = 1:length(seqALPHA)
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i));
end

