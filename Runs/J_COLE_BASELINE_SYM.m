clc
clear
warning off

% cd ../
 
load('VAP31_STEADY_VISCOUS_RELAXED_WING.mat');
vecCD = [OUTP.vecCD]';
vecCLv = [OUTP.vecCLv]';
weight = 13344.6648; % 3000 lb in N
ld = vecCLv./vecCD;
drag = weight./ld;
vecVEHVINF = [OUTP.vecVINF]';

% Finding collective
load('VAP31_CRUISE_PROP_SWEEP_STEADY_FIXED.mat');
temp = [OUTP.vecCT];
vecCT = temp(end,:)';
vecROTORRPM = OUTP(1).vecROTORRPM;
vecROTDIAM = OUTP(1).vecROTDIAM;
thrust = vecCT.*(((vecROTORRPM/60).^2).*((vecROTDIAM).^4))*1.225;
F = scatteredInterpolant([OUTP.vecVINF]', thrust, [OUTP.vecCOLLECTIVE]', 'linear','none');
vecCOLLECTIVE = F(vecVEHVINF, drag./2);


% Not running all the points
% idx = [3 5 6 7];
idx = [5 6 7];
seqALPHA = seqALPHA(idx);
vecVEHVINF = vecVEHVINF(idx);
vecCOLLECTIVE = vecCOLLECTIVE(idx);


clearvars -except seqALPHA vecCOLLECTIVE vecVEHVINF
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    VAP_IN.valSTARTFORCES = 115;
    VAP_IN.valMAXTIME = 120
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end

save('VAP31_STEADY_VISCOUS_FIXED_SYM.mat')
