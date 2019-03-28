clc
clear
warning off
cores = 7;
parpool(cores,'IdleTimeout',800)
filename = 'inputs/WIPP_FINAL.vap';

% Timestep Size Convergence
vecAZNUM = [20,40,60,80,120,160,200];
parfor i = 1:length(vecAZNUM)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 5;
    rps = 4705/60;
    VAP_IN.valDELTIME = 1./(vecAZNUM(i).*rps);
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end

save('TimestepSizeStudy')