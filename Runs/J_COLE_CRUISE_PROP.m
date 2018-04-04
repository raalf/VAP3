clc
clear
warning off

% RPM is 2250
filename = 'inputs/J_COLE_X57_CRUISE_PROP_20_SECTIONS.vap';

% Sweeping collective and speed
vecCOLLECTIVE = [-10:1:10]';
vecVINF = [40:5:90]';
% vecCOLLECTIVE = 0;
% vecVINF = 60;

len = length(vecCOLLECTIVE);
vecCOLLECTIVE = repmat(vecCOLLECTIVE, length(vecVINF), 1);
vecVINF = reshape(repmat(vecVINF, 1, len)', [], 1);

% Timesteps per RPM
timestep_rpm = ceil(1/((2250/60)*0.0005));

for i = 1:length(vecCOLLECTIVE)
    VAP_IN = [];
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
    VAP_IN.vecVEHVINF = vecVINF(i);
    VAP_IN.valSTARTFORCES = 200;
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end

save('J_COLE_CRUISE_PROP_J_CT_SWEEP.mat')