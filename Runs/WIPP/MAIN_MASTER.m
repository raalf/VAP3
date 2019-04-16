clear,clc
warning off
% Setup parpool
% cores = 48;
% parpool(cores,'IdleTimeout',800)


% List of prioity 1 & 2 cases
seqALPHA_PRIORITY1 = [0 5 7 15 17];
seqALPHA_PRIORITY2 = [-7 -5 0 2 5 7 9 11 13 15 17];

seqCT_PRIORITY1 = [0 0.4];
seqCT_PRIORITY2 = [0 0.04 0.2 0.4];

seqMACH_PRIORITY1 = 0.08;
seqMACH_PRIORITY2 = [0.08 0.11];

% Calculate resulting rpm freestream and rpm
seqVEL_PRIORITY1 = seqMACH_PRIORITY1*343;
seqVEL_PRIORITY2 = seqMACH_PRIORITY2*343;

rps_setting = [98.42,78.41,57.6;134.33,107.33,79.33];


for i = 1:length(seqCT_PRIORITY2)
    for j = 1:length(seqMACH_PRIORITY2)
        vel = seqVEL_PRIORITY2(j);
        rps = rps_setting(j,i);
        parfor k = 1:length(seqALPHA_PRIORITY2)
            filename = '/inputs/WIPP_FINAL.vap';
            valAZNUM = 80;
            VAP_IN = [];
            VAP_IN.valMAXTIME = 280;
            VAP_IN.valSTARTFORCES = 200;
            VAP_IN.vecCOLLECTIVE = 6;
            
            VAP_IN.vecVEHALPHA = seqALPHA_PRIORITY2(k);
            VAP_IN.vecVEHVINF = vel;
            VAP_IN.vecROTORRPM = rps*60;
            VAP_IN.valDELTIME = 1./(valAZNUM.*rps);
            OUTP = fcnVAP_MAIN(filename, VAP_IN);
        end
    end
end