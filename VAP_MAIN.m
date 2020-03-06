clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
filename = 'inputs/CREATeV.vap';
% filename = 'inputs/TMotor.vap';
% alpha_seq = -12:2:8;

% for i = 1:length(alpha_seq)
    
%     VAP_IN.vecVEHALPHA = alpha_seq(i);
VAP_IN = [];
[OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(filename, VAP_IN);

% end

