clc
clear

filename = 'inputs/XMLtest.xml';

inp = fcnXML2STRUCT(filename);
VAP = inp.VAP;

clear inp filename

%% Settings
if strcmpi(VAP.settings.flagRELAX.Text, 'true') flagRELAX = 1; else flagRELAX = 0; end
if strcmpi(VAP.settings.flagSTEADY.Text, 'true') flagSTEADY = 1; else flagSTEADY = 0; end
if strcmpi(VAP.settings.flagTRI.Text, 'true') flagTRI = 1; else flagTRI = 0; end

valMAXTIME = int32(str2double(VAP.settings.valMAXTIME.Text));
valMINTIME = int32(str2double(VAP.settings.valMINTIME.Text));
valDELTIME = str2double(VAP.settings.valDELTIME.Text);
valDELTAE = str2double(VAP.settings.valDELTAE.Text);

%% Conditions
valDENSITY = str2double(VAP.conditions.valDENSITY.Text);
valKINV = str2double(VAP.conditions.valKINV.Text);

%% Vehicles
valVEHICLES = max(size(VAP.vehicle));

matVEHORIG = nan(valVEHICLES,3);
vecVEHVINF = nan(valVEHICLES,1);
vecVEHALPHA = nan(valVEHICLES,1);
vecVEHBETA = nan(valVEHICLES,1);
vecVEHROLL = nan(valVEHICLES,1);
vecVEHPITCH = nan(valVEHICLES,1);
vecVEHYAW = nan(valVEHICLES,1);

vecWINGS = nan(valVEHICLES,1);

for i = 1:valVEHICLES
    
    veh = VAP.vehicle{1,i};
    matVEHORIG(i,:) = [str2double(veh.x.Text) str2double(veh.y.Text) str2double(veh.z.Text)];
    vecVEHVINF(i,1) = str2double(veh.vinf.Text);
    vecVEHALPHA(i,1) = str2double(veh.alpha.Text);
    vecVEHBETA(i,1) = str2double(veh.beta.Text);
    vecVEHROLL(i,1) = str2double(veh.roll.Text);
    vecVEHPITCH(i,1) = str2double(veh.pitch.Text);
    vecVEHYAW(i,1) = str2double(veh.yaw.Text);
    
    vecWINGS(i,1) = max(size(veh.wing));
    
    for j = 1:vecWINGS
    
        
        
    end
    
end
