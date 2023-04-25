clc
clear
warning off

cd('../')

addpath('Flight Dynamics')

filename = 'inputs/Discus2c.vap';

trim_iter = 1;

VAP_IN = [];
TRIM = [];

% Initialize variables and read in geometry
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN);
% load('C:\Users\Michael\Desktop\Optimum_SS.mat')
FLAG.OPT = 1; % Flag to load existing optimization data (1), parametric sweep data (2) or use input file data (0)
% FLAG.STRUCTURE = 0;
COND.valMAXTRIMITER = 50;

% Parametric sweep data to load (assumes constant properties across the
% span)
% INPU.matEIx_param = 1500000;
% INPU.matGJt_param = 100000;
% INPU.vecEA_param = 0.4;
% INPU.vecCG_param = 0.4;
% INPU.matEIx_param = 100000;
% INPU.matGJt_param = 250000;
% INPU.vecEA_param = 0.4;
% INPU.vecCG_param = 0.4;

% load('C:\Users\Michael\Desktop\Optimum_SS.mat')
% OUTP.matDEF = 0*OUTP.matDEF;
% OUTP.matDEFGLOB = 0*OUTP.matDEFGLOB;
% OUTP.matTWIST = 0*OUTP.matTWIST;
% OUTP.matTWISTGLOB = 0*OUTP.matTWISTGLOB;

% Load in optimization case if wanting to run locally
opthistory = importdata('optimization_test.txt');

des = 2;
INPU.matEIx(:,1) = opthistory(des,3:21);
INPU.matGJt(:,1) = opthistory(des,22:40);
INPU.vecEA(:,1) = opthistory(des,41:59);
INPU.vecCG(:,1) = opthistory(des,60:78);

% Initialize DVE parameters
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP);

% Fix vehicle CG to a specific x-location (comment 47-65 out if not fixing
% CG location)
CGX_loc = 0.4;
Xnew = (CGX_loc*COND.vecVEHWEIGHT/9.81 - VEHI.vecWINGMASS(1)*VEHI.vecWINGCG(1,1))/sum(VEHI.vecFUSEMASS,1);
[INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle

VEHI.vecFUSELOC = [Xnew - VEHI.vecFUSEL/2,0,0];
VEHI.vecFUSEMASS = sum(VEHI.vecFUSEMASS,1); % Turn fuse mass from distributed to point mass
VEHI.vecFUSELM = VEHI.vecFUSEMASS/VEHI.vecFUSEL;
VEHI.valFUSEDX = VEHI.vecFUSEL/(VEHI.valNFELE-1);
tempdx = [[0; cumsum(repmat(VEHI.valFUSEDX,VEHI.valNFELE-1,1))],zeros(VEHI.valNFELE,1),zeros(VEHI.valNFELE,1);];
VEHI.vecFUSEBEAM = repmat(VEHI.vecFUSELOC,VEHI.valNFELE,1) + tempdx;

VEHI.vecFUSEMASSLOC = (VEHI.vecFUSEBEAM(2:end,:) + VEHI.vecFUSEBEAM(1:end-1,:))./2;
tempFUSELM = interp1(VEHI.vecFUSEBEAM(:,1),repmat(VEHI.vecFUSELM,VEHI.valNFELE,1),VEHI.vecFUSEMASSLOC(:,1));
VEHI.vecFUSEMASS = tempFUSELM.*VEHI.valFUSEDX;
VEHI.vecFUSECG = sum(VEHI.vecFUSEMASS.*(VEHI.vecFUSEMASSLOC-INPU.matVEHORIG),1)./(sum(VEHI.vecFUSEMASS,1)); % Fuse CG location relative to vehicle origin

VEHI.vecFUSEBEAM = fcnGLOBSTAR(VEHI.vecFUSEBEAM, zeros(VEHI.valNFELE,1), repmat(deg2rad(-COND.vecVEHPITCH),VEHI.valNFELE,1), zeros(VEHI.valNFELE,1));
VEHI.vecFUSEMASSLOC = fcnGLOBSTAR(VEHI.vecFUSEMASSLOC, zeros(VEHI.valNFELE-1,1), repmat(deg2rad(-COND.vecVEHPITCH),VEHI.valNFELE-1,1), zeros(VEHI.valNFELE-1,1));
[INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle

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
FLAG.VISCOUS = 1; % Turn viscous calcs on (1) or off (0)
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

[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);

fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)

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
fprintf('\nVehicle trimmed. AoA = %.3f deg., Elev. Angle = %.3f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)
fprintf('\nDeflection = %.3f m, Twist = %.4f deg.\n\n',OUTP.matDEFGLOB(end,end),rad2deg(OUTP.matTWISTGLOB(end,end)))
end

% load('Discus_Optimum_Trim.mat')
FLAG.STIFFWING = 0;
FLAG.GIF = 0;
%% Perform full flight-dynamic simulation on trimmed/deformed aircraft
SURF.matBEAMACC = [];
COND.valGUSTAMP = 1;
COND.valGUSTSTART = 40;
OUTP.valVS = COND.vecVEHVINF*sind(COND.vecVEHFPA); % Steady-state sink rate

SURF.matB = 0*[max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];

COND.valMAXTIME = ceil((10*COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART); % Max time is 5 gust lengths
% COND.valMAXTIME = 750;
COND.valSTIFFSTEPS = 15; % Don't compute structure deflections for 15 timesteps
COND.valSTARTFORCES = 1; % Start force calcs at timestep 1
FLAG.FLIGHTDYN = 1; % Flight dynamics model on (1) or off (0)
FLAG.STATICAERO = 0;
FLAG.STEADY = 0; % Steady (1) or unsteady (2) aerodynamics
FLAG.RELAX = 0; % Don't set this to 1 unless you want to have a bad time
FLAG.GUSTMODE = 4; % Sine (1), 1-cosine (2), sharp-edge (3), von Karman (4) [MUST LOAD GUST PROFILE IF YOU DON'T WANT IT TO CHANGE EVERY RUN] or no gust (0)
FLAG.SAVETIMESTEP = 0;

VEHI.vecVEHDYN(1:COND.valSTIFFSTEPS,4) = deg2rad(COND.vecVEHPITCH); % Initial conditions for flight-dynamic model

[VEHI.matVEHUVW] = fcnGLOBSTAR(VEHI.matGLOBUVW, 0, pi+deg2rad(COND.vecVEHPITCH), 0);

COND.start_loc = repmat([-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF,0,0],size(SURF.matCENTER,1),1, size(SURF.matCENTER,3)); % Location (in meters) in global frame where gust starts

%% Load continuous gust profile
load('vKarman_Gust.mat');
SURF.vk_gust = [cgust; 0];
temp = 0:COND.valDELTIME*COND.vecVEHVINF:(length(cgust)*COND.vecVEHVINF*COND.valDELTIME);
SURF.matGUSTFIELD = [-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF - temp',zeros(length(temp),2)];
SURF.vk_gust = repmat(SURF.vk_gust,1,size(SURF.matCENTER,1));
SURF.vk_gust = SURF.vk_gust';
SURF.vk_gust = reshape(SURF.vk_gust,size(SURF.matCENTER,1),1,length(temp));

%% Timestep
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);
