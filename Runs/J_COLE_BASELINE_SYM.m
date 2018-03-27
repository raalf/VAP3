clc
clear
warning off

% cd ../
 
% Approximating L/D based on OVERFLOW
load('OVERFLOW_TRANSITION.mat');
ld = OVERFLOW_TRANSITION(:,2)./OVERFLOW_TRANSITION(:,1);
weight = 13344.6648; % 3000 lb in N
drag = weight./ld;

% Finding collective
load('VAP315_STEADY_VISCOUS_FIXED_CRUISE_PROP_20RadialSections')
temp = [OUTP(:).vecCT];
vecCT = temp(end,14:22)';
vecROTORRPM = 1719;
vecROTDIAM = 1.524;
thrust = vecCT.*(((vecROTORRPM/60).^2).*((vecROTDIAM).^4))*1.225;
vecCOLLECTIVE = interp1(thrust, vecCOLLECTIVE(14:22)', drag./2)+20;

% Finding alpha
load('VAP315_STEADY_VISCOUS_RELAXED_WING.mat','OUTP', 'seqALPHA')
vecCD = [OUTP.vecCD]';
vecCLv = [OUTP.vecCLv]';
cl = interp1(vecCD, vecCLv, OVERFLOW_TRANSITION(:,1));
seqALPHA = interp1(vecCLv, seqALPHA, cl);

clearvars -except seqALPHA vecCOLLECTIVE
filename = 'inputs/J_COLE_BASELINE_SYM.vap';
parfor i = 1:length(vecCOLLECTIVE)
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), vecCOLLECTIVE(i));
end

save('VAP31_STEADY_VISCOUS_FIXED_SYM.mat')
