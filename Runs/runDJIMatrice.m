clc
clear
%% Sort out folder and paths
folder = pwd;
if strcmp(folder(end-4:end),'\Runs')
    addpath(folder)
    cd(folder(1:end-4))
end

%% Run fcnVAP over a series of advance ratios and angles
angles = [-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90];
rpm = [2000 3000 4000 5000 6000];
J = 0:0.1:1;

filename = 'inputs/Matrice_210_RTK_Rotor.vap';

OUTP = fcnVAP_MAIN(filename, []);

parfor i = 1:length(angles)
    AOA_OI = angles(i);
    
    for j = 1:length(J)
        adv_r = J(j);
        
        for k = 1:length(rpm)
            VAP_IN = [];
            VAP_IN.valMAXTIME = 160;
            VAP_IN.valSTARTFORCES = 140;
            VAP_IN.vecROTORRPM = rpm(k);
            VAP_IN.vecVEHALPHA = AOA_OI;
            VAP_IN.vecVEHVINF = adv_r*(rpm(k)/60)*0.4572;
            VAP_IN.valDELTIME = 1/((rpm(k)/60)*20); %20 timesteps per rev
            
            OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
        end
    end
end