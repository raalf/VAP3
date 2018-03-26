clc
clear

seqALPHA = [2:1:14];

seqALPHA = 10

filename = 'inputs/CREATeV.vap';
for i = 1:length(seqALPHA)
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), 0);
end