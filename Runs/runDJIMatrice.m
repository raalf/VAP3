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
J = 0.1:0.1:1;

filename = 'inputs/Matrice_210_RTK_Rotor.vap';
CT = nan(length(J),length(angles),length(rpm));
CP = nan(length(J),length(angles),length(rpm));
CQ = nan(length(J),length(angles),length(rpm));
CMx = nan(length(J),length(angles),length(rpm));
CMy = nan(length(J),length(angles),length(rpm));
CNx = nan(length(J),length(angles),length(rpm));
CNy = nan(length(J),length(angles),length(rpm));

leni = length(angles);
lenj = length(J);
lenk = length(rpm);
for i = 1:leni
    AOA_OI = angles(i);
    
    for j = 1:lenj
        adv_r = J(j);
        
        for k = 1:lenk
            VAP_IN = [];
            VAP_IN.RELAX = false;
            VAP_IN.valMAXTIME = 160;
            VAP_IN.valSTARTFORCES = 140;
            VAP_IN.vecROTORRPM = rpm(k);
            VAP_IN.vecVEHALPHA = AOA_OI;
            VAP_IN.vecVEHVINF = adv_r*(rpm(k)/60)*0.4572;
            VAP_IN.valDELTIME = 1/((rpm(k)/60)*20); %20 timesteps per rev
            
            OUTP = fcnVAP_MAIN(filename, VAP_IN);
            CT(j,i,k) = OUTP.vecCT_AVG;
            CP(j,i,k) = OUTP.vecCP_AVG;
            CQ(j,i,k) = OUTP.vecCP_AVG/(2*pi*rpm(k)/60);
            CMx(j,i,k) = OUTP.vecCMx_AVG;
            CMy(j,i,k) = OUTP.vecCMy_AVG;
            CNx(j,i,k) = OUTP.vecCFx_AVG;
            CNy(j,i,k) = OUTP.vecCFy_AVG;
        end
    end
end

[Angle, J, RPM] = meshgrid(angles,J,rpm);

save('DJI_Matrice_Data')