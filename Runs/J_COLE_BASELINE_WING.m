clc
clear
warning off

% load('VAP31_STEADY_VISCOUS_RELAXED_WING.mat');
% vecCD = [OUTP.vecCD]';
% vecCLv = [OUTP.vecCLv]';
% clearvars -except vecCD vecCLv

seqALPHA = [2:2:14];

filename = 'inputs/J_COLE_BASELINE_WING.vap';
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.RELAX = 0;
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end
