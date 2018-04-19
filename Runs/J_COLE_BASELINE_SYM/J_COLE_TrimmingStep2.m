clear
clc

% Import 1st Iteration Data
load('VAP32_WING+PROP_forCDo_TS160_fixed_last22Timesteps.mat')
OUTP1 = OUTP;
clear OUTP;

CL1 = nanmean([OUTP1.vecCL],1);
CT1 = nanmean([OUTP1.vecCT],1);

% delta of targeted CL, CT and calculated CL, CT from 1st iteration
dCL1 = CL - CL1;
dCT1 = CT - CT1;

% interpolate alpha to maintain steady level flight at VINF 
% using wing only data
% New sets of AOA input for 2nd iteration in order to hit the targeted CL
seqALPHA2 = interp1([WING.OUTP.vecCLv],[WING.OUTP.vecVEHALPHA],CL + dCL1);
% New sets of collective pitch input for 2nd iteration in order to hit the targeted CT
vecCOLLECTIVE2 = F(vecVEHVINF, CT + dCT1);



% Running
warning off
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE2)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA2(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    VAP_IN.valSTARTFORCES = 138;
    VAP_IN.valMAXTIME = 160;
    
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);

    fprintf('finished AOA=%.1f\n',seqALPHA(i));
end

save('VAP32_WPforCDo_TS160_22_fixed_iter2.mat')




