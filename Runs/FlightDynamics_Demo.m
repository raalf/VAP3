clc
clear
warning off

cd('../')

addpath('Flight Dynamics')

filename = 'inputs/Flight_Dynamics.vap';

trim_iter = 1;

VAP_IN = [];
TRIM = [];

% Initialize variables and read in geometry
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN);
FLAG.OPT = 0; % Flag to load existing optimization data (1), parametric sweep data (2) or use input file data (0)
COND.valMAXTRIMITER = 50; % Max trim iterations

% Initialize DVE parameters
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP);

COND.valSTIFFSTEPS = inf;

% Required CL for L = W for given speed
VEHI.vecVEHPITCH = COND.vecVEHALPHA; % Set vehicle pitch angle equal to AoA for initial time stepping

FLAG.FLIGHTDYN = 0;
OUTP.dyn_iter = 0;
FLAG.PLOT = 0;

% Initialize trim tolerance variables
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

[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC); % Trim iteration function

fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)

FLAG.STIFFWING = 0; % Turn on flex wing calcs

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
COND.valGUSTSTART = 40; % Time step gust starts at
OUTP.valVS = COND.vecVEHVINF*sind(COND.vecVEHFPA); % Steady-state sink rate

SURF.matB = 0*[max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4]; % Structural damping; set to zero for now

valGUSTT = 1; % Number of gust periods for simulation (will only run one gust plus valGUSTT - 1 gust lengths of no gust)
COND.valMAXTIME = ceil((valGUSTT*COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART); % Max time is valGUSTT gust periods
% COND.valMAXTIME = 750;
COND.valSTIFFSTEPS = 15; % Don't compute structure deflections for 15 timesteps
COND.valSTARTFORCES = 1; % Start force calcs at timestep 1
FLAG.FLIGHTDYN = 1; % Flight dynamics model on (1) or off (0)
FLAG.STATICAERO = 0;
FLAG.STEADY = 0; % Steady (1) or unsteady (0) aerodynamics
FLAG.RELAX = 0; % Don't set this to 1 unless you want to have a bad time
FLAG.GUSTMODE = 1; % Sine (1), 1-cosine (2), sharp-edge (3), von Karman (4) [MUST LOAD GUST PROFILE IF YOU DON'T WANT IT TO CHANGE EVERY RUN] or no gust (0)
FLAG.SAVETIMESTEP = 0;

VEHI.vecVEHDYN(1:COND.valSTIFFSTEPS,4) = deg2rad(COND.vecVEHPITCH); % Use trimmed pitch angle for initial conditions for flight-dynamic model

[VEHI.matVEHUVW] = fcnGLOBSTAR(VEHI.matGLOBUVW, 0, pi+deg2rad(COND.vecVEHPITCH), 0); % Compute body frame U, V, W velocities

COND.start_loc = repmat([-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF,0,0],size(SURF.matCENTER,1),1, size(SURF.matCENTER,3)); % Location (in meters) in global frame where gust starts

%% Load continuous gust profile (comment out if using discrete gusts)
load('vKarman_Gust.mat');
SURF.vk_gust = [cgust; 0];
temp = 0:COND.valDELTIME*COND.vecVEHVINF:(length(cgust)*COND.vecVEHVINF*COND.valDELTIME);
SURF.matGUSTFIELD = [-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF - temp',zeros(length(temp),2)];
SURF.vk_gust = repmat(SURF.vk_gust,1,size(SURF.matCENTER,1));
SURF.vk_gust = SURF.vk_gust';
SURF.vk_gust = reshape(SURF.vk_gust,size(SURF.matCENTER,1),1,length(temp));

%% Timestep
[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);
