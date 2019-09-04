% clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
% filename = 'inputs/TMotor.vap';
filename = 'inputs/Goland_Wing.vap';
% alpha_seq = -5:10;

% for i = 1:length(alpha_seq)
    
    VAP_IN = [];
    [OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(filename, VAP_IN);
%     CL(i) = OUTP.vecCL(end);
%     CDi(i) = OUTP.vecCDI(end);
    
% end
save('Sharp_Edge.mat')