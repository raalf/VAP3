% clc
clear

cd('C:\Users\travi\OneDrive\Desktop\GIT\VAP3\Runs\Winglet Optimization')

if exist('aux_files','file') ~= 7
    mkdir('aux_files');
end

% relax, start force, aux_files

z = [0.5 ...
    0 7.45 0.53 0.2 1 ...
    0 7.5 0.55 0.1 4 ...
    0.2 7.46 0.6 0.2 1 ...
    0.3 7.5 0.9 0.1 -2];

%% Setting up files
vap_filename = ['aux_files/', regexprep(tempname('/'), {'/', '\'}, ''), '.vap'];
copyfile('Standard_Cirrus.vap', vap_filename);

geom = [0 0 0 0.92 0; 0 4.5 0.315 0.644 0.3; 0 7.5 0.525 0.35 1.8];

geom(3,2) = 7.5 - z(1);
geom(3,3) = geom(2,3) + (7.5 - (z(1)) - 4.5).*((geom(3,3) - geom(2,3))./(geom(3,2) - geom(2,2)));
geom(3,4) = geom(2,4) + (7.5 - (z(1)) - 4.5).*((geom(3,4) - geom(2,4))./(geom(3,2) - geom(2,2)));

pan(1).geom = geom;
pan(2).geom = [geom(end,:); z(2:6); z(7:11)];
pan(3).geom = [geom(end,:); z(12:16); z(17:end)];

vap3_inputmod_wing(vap_filename, pan)

cd ./../../
seqALPHA = 10;
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.RELAX = 0;
    VAP_IN.valSTARTFORCES = 40;
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    WING_SWEEP(i) = fcnVAP_MAIN(['Runs/Winglet Optimization/', vap_filename], VAP_IN);
end
cd('Runs/Winglet Optimization/');
delete(vap_filename)

%% Analysis
% Calculate CL at VINF and S&L flight
S    = WING_SWEEP(1).valAREA; % ref. platform area
CL   = [WING_SWEEP.vecCLv_AVG]';

%% Output




