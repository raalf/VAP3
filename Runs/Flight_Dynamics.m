clc
clear
warning off

addpath('Flight Dynamics')

filename = 'inputs/Discus2c.vap';

trim_iter = 1;

VAP_IN = [];
TRIM = [];

% Initialize variables and read in geometry
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN);

FLAG.OPT = 2;
COND.valMAXTRIMITER = 50;

INPU.matEIx_param = [100000; 100000; 1500000];
INPU.matGJt_param = [500000; 500000; 100000];
INPU.vecEA_param = [0.65; 0.65; 0.65];
INPU.vecCG_param = [0.25; 0.25; 0.25];

% opthistory = importdata('G:\My Drive\PhD\Optimization\opthistory_cosine.txt');

% des = 15;
% INPU.matEIx(:,1) = opthistory(des,2:20);
% INPU.matGJt(:,1) = opthistory(des,21:39);
% INPU.vecEA(:,1) = opthistory(des,40:58);
% INPU.vecCG(:,1) = opthistory(des,59:77);

[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP);

COND.valSTIFFSTEPS = inf;

% Required CL for L = W for given speed
COND.CLtrim = 2*COND.vecVEHWEIGHT/(COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
VEHI.vecVEHPITCH = COND.vecVEHALPHA;

FLAG.FLIGHTDYN = 0;
OUTP.dyn_iter = 0;
FLAG.PLOT = 0;
tol1 = 100;
tol2 = 100;
tol3 = 100;
tol = 100;

% Initial trim iteration
% Find trimmed conditions for rigid aircraft. These loads are used to
% compute the static aeroelastic deflections

FLAG.TRIM = 0;
FLAG.VISCOUS = 0;
FLAG.GUSTMODE = 0;

% Increase number of stiff steps (steps before computing structural
% deflection) to happen beyond valMAXTIME. Keeps solution to be for a rigid
% aircraft
COND.valSTIFFSTEPS = COND.valMAXTIME + 1;

COND.valSTARTFORCES = COND.valMAXTIME; % Only compute forces on last timestep

SURF.vecELEVANGLE = 0;

SURF.center_dist = cumsum((2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)))-SURF.vecDVEHVSPN(1);

% Initial trim iteration loop for rigid aircraft
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 0);
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

SURF.center_dist = cumsum((2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)))-SURF.vecDVEHVSPN(1);
VEHI.vecPROPLOC_START = VEHI.vecPROPLOC;

% fun = @(x)fcnTRIMOPT(x, FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
% options = optimset('TolFun',1e-2);
% nvars = 2;
% x = fminsearch(fun,[5;0],options);

[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

% tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)

% save('HALE_Validation_10ms_Rigid.mat')
FLAG.STIFFWING = 0;

OUTP.aero_iter = 1;
if FLAG.STIFFWING == 0
FLAG.TRIM = 0;
FLAG.STATICAERO = 1;
OUTP.aero_iter = 1;
    
%% Trim iteration loop 
% Update aeroelastic deflection based on trim conditions and then update
% trim conditions
while max(abs(tol)) > 0.01
    COND.valSTIFFSTEPS = COND.valMAXTIME - 1;
    
    COND.valSTARTFORCES = COND.valMAXTIME - 2; % Compute forces for last two timesteps
%     COND.valSTARTFORCES = COND.valMAXTIME; % Compute forces for last two timesteps
    
%     COND.valSTARTFORCES = 1;
    
    tol_aero = 100;
    tol_aero2 = 100;
    tol_aero1 = 100;
    SURF.matB = [max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];
    
    %% Aeroelastic convergence loop
    % Compute static aeroelastic configuration based on trimmed aircraft loads
    while max(abs(tol_aero)) > 1e-3
        
        [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 0);

        [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
        
        % Compute tolerance for solution convergence. Tip deflection and
        % tip twist are compared to their previous iterations solution
        temp_def(OUTP.aero_iter,1) = OUTP.matDEFGLOB(end,end);
        temp_twist(OUTP.aero_iter,1) = OUTP.matTWISTGLOB(end,end);
        
        if OUTP.aero_iter > 1
            tol_aero1 = (temp_def(OUTP.aero_iter,1)-temp_def(OUTP.aero_iter-1,1))/temp_def(OUTP.aero_iter-1,1);
            tol_aero2 = (temp_twist(OUTP.aero_iter,1)-temp_twist(OUTP.aero_iter-1,1))/temp_twist(OUTP.aero_iter-1,1);
        end
        
        tol_aero = [tol_aero1; tol_aero2];
        
        OUTP.aero_iter = OUTP.aero_iter + 1;
    end    
    
    %% Trim new deformed configuration
    COND.valSTIFFSTEPS = COND.valMAXTIME + 1;
    COND.valSTARTFORCES = COND.valMAXTIME;
    
    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 0);
    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    % Compute tolerance for trim solution convergence. Trim angle of
    % attack and wing tip deflection are compared to their previous
    % iteration
    trim_def(trim_iter,1) = OUTP.matDEFGLOB(end,end);
    trim_twist(trim_iter,1) = OUTP.matTWISTGLOB(end,end);
    trim_alpha(trim_iter,1) = COND.vecVEHALPHA;

    if trim_iter > 1
        tol1 = (trim_alpha(trim_iter,1)-trim_alpha(trim_iter-1,1))/trim_alpha(trim_iter-1,1);
        tol2 = (trim_def(trim_iter,1)-trim_def(trim_iter-1,1))/trim_def(trim_iter-1,1);
        tol3 = (trim_twist(trim_iter,1)-trim_twist(trim_iter-1,1))/trim_twist(trim_iter-1,1);
    end
    
    tol = [tol1; tol2; tol3];

    trim_iter = trim_iter + 1;
end

OUTP.aero_iter = 0;
tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)
end
% save('Discus_Rigid_Trim.mat')
save('temp_Flex_Trim.mat')
%% Perform full flight-dynamic simulation on trimmed/deformed aircraft
% load('HALE_Flex_Trim.mat')
SURF.matBEAMACC = [];
COND.valGUSTAMP = 1;
COND.valGUSTL = 50;
COND.valGUSTSTART = 40;

SURF.matB = [max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];

% valTBOOM = SURF.matVLST(SURF.matDVE(SURF.idxTAIL(1),1),1) - SURF.matVLST(SURF.matDVE(SURF.idxFLEX(1),1),1);
COND.valMAXTIME = ceil((COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART);
% COND.valMAXTIME = 465;
COND.valSTIFFSTEPS = 15;
COND.valSTARTFORCES = 1;
FLAG.FLIGHTDYN = 1;
FLAG.STATICAERO = 0;
FLAG.STEADY = 0;
FLAG.RELAX = 0;
FLAG.GUSTMODE = 2;
FLAG.SAVETIMESTEP = 0;

VEHI.vecVEHDYN(1:COND.valSTIFFSTEPS,4) = deg2rad(COND.vecVEHPITCH);

[VEHI.matVEHUVW] = fcnGLOBSTAR(VEHI.matGLOBUVW, 0, pi+deg2rad(COND.vecVEHPITCH), 0);

COND.start_loc = repmat([-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF,0,0],size(SURF.matCENTER,1),1, size(SURF.matCENTER,3)); % Location (in meters) in global frame where gust starts

% OUTP.vecFUSEREACTION = sum(OUTP.vecBEAMFORCE(2:end).*INPU.valDY,1);

[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);

save('C:\Users\Michael\Desktop\Grid_Density4.mat');

temp_gain = OUTP.vecZE_old(end);

clearvars -except temp_gain
load('temp_Flex_Trim.mat')
SURF.matBEAMACC = [];
COND.valGUSTAMP = 1;
COND.valGUSTL = 50;
COND.valGUSTSTART = 40;

SURF.matB = [max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];

COND.valMAXTIME = ceil((COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART);
COND.valSTIFFSTEPS = 15;
COND.valSTARTFORCES = 1;
FLAG.FLIGHTDYN = 1;
FLAG.STATICAERO = 0;
FLAG.STEADY = 0;
FLAG.RELAX = 0;
FLAG.GUSTMODE = 0;
FLAG.SAVETIMESTEP = 0;
FLAG.STIFFWING = 1;

VEHI.vecVEHDYN(1:COND.valSTIFFSTEPS,4) = deg2rad(COND.vecVEHPITCH);

[VEHI.matVEHUVW] = fcnGLOBSTAR(VEHI.matGLOBUVW, 0, pi+deg2rad(COND.vecVEHPITCH), 0);

COND.start_loc = repmat([-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF,0,0],size(SURF.matCENTER,1),1, size(SURF.matCENTER,3)); % Location (in meters) in global frame where gust starts


[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);

gain = temp_gain - OUTP.vecZE_old(end);

delete temp_Flex_Trim.mat;
