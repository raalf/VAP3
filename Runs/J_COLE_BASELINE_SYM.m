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
load('VAP315_STEADY_VISCOUS_FIXED_CRUISE_PROP_20RadialSections')
temp = [OUTP(:).vecCT];
vecCT = temp(end,14:22)';
vecROTORRPM = 1719;
vecROTDIAM = 1.524;
thrust = vecCT.*(((vecROTORRPM/60).^2).*((vecROTDIAM).^4))*1.225;
vecCOLLECTIVE = interp1(thrust, vecCOLLECTIVE(14:22)', drag./2) + 20;

clearvars -except seqALPHA vecCOLLECTIVE vecVEHVINF
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
for i = 1:length(vecCOLLECTIVE)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
    VAP_IN.vecVEHVINF = vecVEHVINF(i);
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), vecCOLLECTIVE(i));
end

% save('VAP31_STEADY_VISCOUS_FIXED_SYM.mat')
