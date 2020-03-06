clc
clear

seqALPHA = [-4:1:8];

% seqALPHA = 10

filename = 'inputs/CREATeV.vap';
% for i = 1:length(seqALPHA)
    OUTP = fcnVAP_MAIN(filename, []);
% end