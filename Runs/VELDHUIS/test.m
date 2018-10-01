clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/VELDHUIS/aux_files')
addpath('./../../Runs/VELDHUIS/')

% NACA 64_2 A015
% Prop is 0.2018 m in front of wing
% Prop is 4 bladed, diameter 0.236 m, 3/4 R pitch at 25 degrees

prop_y = [0.24:0.08:0.64];

OUT = [];

%%
parfor i = 1:size(prop_y,2)
    rotor = [];
    
    % Constant propeller values
    rotor.hub = [-0.2018 prop_y(i) 0];
    rotor.dir = 1;
    rotor.collective = 0;
    rotor.axis = [-1 0 0];
    rotor.m = 1;
    rotor.blades = 4;
    rotor.airfoil = 'MH-117';
    rotor.diam = 0.236;
    airfoil_data = load('airfoils/MH-117.mat');
    
    % Conditions
    rho = 1.225;
    Re = 3e+5;
    kinv = 1.45e-5;
    cmac = 0.24;
    vinf = (Re*kinv)/cmac;
    J = 0.92;
    n = vinf/(J*rotor.diam);
    rotor.rpm = n*60;
    T_C = 0.025; % (T/(rho.*V.^2.*D.^2))
    alpha = 4.2;
    thrust = T_C*rho*(n^2)*(rotor.diam^2);
    
    
    %% QMIL
    temp_name = regexprep(tempname('/'), {'/', '\'}, '');
    qmil_path = fcnQMILCREATE(temp_name, airfoil_data, rotor.blades, thrust, vinf, rotor.rpm, rotor.diam);
    qmil_filename = regexprep(qmil_path, 'aux_files/', '');
    qmil_output_filename = ['output_', qmil_filename];
    if ispc
        qmil_output_path = ['aux_files\', qmil_output_filename];
        exeName = ['aux_files\', qmil_filename, '.exe'];
        copyfile('qmil.exe', exeName);
    else
        qmil_output_path = ['aux_files/', qmil_output_filename];
        exeName = ['aux_files/', qmil_filename, 'ex'];
        copyfile('qmilex', exeName);
    end
    
    prmpt = sprintf('%s %s %s', exeName, qmil_path, qmil_output_path);
    [~,~] = system(prmpt);
    
    delete(exeName);
    delete(qmil_path);
    
    %% Building Propeller & Wing VAP File
    vap_filename = ['aux_files/vap_',temp_name,'.vap'];
    copyfile('VELDHUIS_BASELINE_WING.vap', vap_filename);
    vap3_inputmod_prop(vap_filename, rotor, qmil_output_path);
    
    delete(qmil_output_path);
    
    %%
    cd '../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = alpha;
    VAP_IN.vecVEHVINF = vinf;
    VAP_IN.valMAXTIME = 160;
    VAP_IN.valSTARTFORCES = VAP_IN.valMAXTIME-20;
%                     VAP_IN.valMAXTIME = 5
%                     VAP_IN.valSTARTFORCES = 4
    VAP_IN.valDELTIME = (1/60)/(rotor.rpm/60);
    OUTP = fcnVAP_MAIN(vap_filename, VAP_IN);
    
    OUT(i).OUTP = OUTP;
    
    cd 'Runs/VELDHUIS/'
    
    delete(vap_filename)
end

%%
load('sweep.mat');

hFig1 = figure(1);
clf(1);

hold on
for i = 1:size(prop_y,2)
    scatter(prop_y(i)/0.64, OUT(i).OUTP.vecCL_AVG/OUT(i).OUTP.vecCD_AVG,'sk')
end
box on
grid minor
axis tight

xlabel('Propeller Hub Location (y/b/2)','FontSize',15);
ylabel('C_L/C_D','FontSize',15);












