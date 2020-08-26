clc
clear
warning off

addpath('Flight Dynamics')

% filename = 'inputs/HALE_new.vap';
% 
% trim_iter = 1;
% 
% VAP_IN = [];
% TRIM = [];
% 
% % Initialize variables and read in geometry
% [FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN);
% 
% [FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP);
% 
% COND.valSTIFFSTEPS = inf;
% 
% % Required CL for L = W for given speed
% COND.CLtrim = 2*COND.vecVEHWEIGHT/(COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA);
% 
% FLAG.FLIGHTDYN = 0;
% OUTP.dyn_iter = 0;
% FLAG.PLOT = 0;
% tol1 = 100;
% tol2 = 100;
% tol = 100;
% 
% % Initial trim iteration
% % Find trimmed conditions for rigid aircraft. These loads are used to
% % compute the static aeroelastic deflections
% 
% FLAG.TRIM = 0;
% FLAG.VISCOUS = 0;
% 
% % Increase number of stiff steps (steps before computing structural
% % deflection) to happen beyond valMAXTIME. Keeps solution to be for a rigid
% % aircraft
% COND.valSTIFFSTEPS = COND.valMAXTIME + 1;
% 
% COND.valSTARTFORCES = COND.valMAXTIME; % Only compute forces on last timestep
% 
% SURF.vecELEVANGLE = 0;
% 
% SURF.center_dist = cumsum((2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)))-SURF.vecDVEHVSPN(1);
% 
% % Initial trim iteration loop for rigid aircraft
% [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);
% [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
% 
% SURF.center_dist = cumsum((2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)))-SURF.vecDVEHVSPN(1);
% VEHI.vecPROPLOC_START = VEHI.vecPROPLOC;
% 
% % fun = @(x)fcnTRIMOPT(x, FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
% % options = optimset('TolFun',1e-2);
% % nvars = 2;
% % x = fminsearch(fun,[5;0],options);
% 
% [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
% 
% tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
% fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,tail_angle)
% 
% save('HALE_Validation_10ms_Rigid.mat')
% 
% FLAG.TRIM = 0;
% OUTP.aero_iter = 1;
%     
% %% Trim iteration loop 
% % Update aeroelastic deflection based on trim conditions and then update
% % trim conditions
% while max(abs(tol)) > 0.01
%     COND.valSTIFFSTEPS = COND.valMAXTIME - 1;
%     
%     COND.valSTARTFORCES = COND.valMAXTIME - 2; % Compute forces for last two timesteps
% %     COND.valSTARTFORCES = COND.valMAXTIME; % Compute forces for last two timesteps
%     
% %     COND.valSTARTFORCES = 1;
%     
%     tol_aero = 100;
%     tol_aero2 = 100;
%     tol_aero1 = 100;
%     
%     %% Aeroelastic convergence loop
%     % Compute static aeroelastic configuration based on trimmed aircraft loads
%     while max(abs(tol_aero)) > 1e-3
%         
%         [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
% 
%         [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
%         
%         % Compute tolerance for solution convergence. Tip deflection and
%         % tip twist are compared to their previous iterations solution
%         temp_def(OUTP.aero_iter,1) = OUTP.matDEFGLOB(end,end);
%         temp_twist(OUTP.aero_iter,1) = OUTP.matTWISTGLOB(end,end);
%         
%         if OUTP.aero_iter > 1
%             tol_aero1 = (temp_def(OUTP.aero_iter,1)-temp_def(OUTP.aero_iter-1,1))/temp_def(OUTP.aero_iter-1,1);
%             tol_aero2 = (temp_twist(OUTP.aero_iter,1)-temp_twist(OUTP.aero_iter-1,1))/temp_twist(OUTP.aero_iter-1,1);
%         end
%         
%         tol_aero = [tol_aero1; tol_aero2];
%         
%         OUTP.aero_iter = OUTP.aero_iter + 1;
%     end    
%     
%     %% Trim new deformed configuration
%     COND.valSTIFFSTEPS = COND.valMAXTIME + 1;
%     COND.valSTARTFORCES = COND.valMAXTIME;
%     
%     [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
%     [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnTRIMITER(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
%     
%     % Compute tolerance for trim solution convergence. Trim angle of
%     % attack and wing tip deflection are compared to their previous
%     % iteration
%     trim_def(trim_iter,1) = OUTP.matDEFGLOB(end,end);
%     trim_alpha(trim_iter,1) = COND.vecVEHALPHA;
% 
%     if trim_iter > 1
%         tol1 = (trim_alpha(trim_iter,1)-trim_alpha(trim_iter-1,1))/trim_alpha(trim_iter-1,1);
%         tol2 = (trim_def(trim_iter,1)-trim_def(trim_iter-1,1))/trim_def(trim_iter-1,1);
%     end
%     
%     tol = [tol1; tol2];
% 
%     trim_iter = trim_iter + 1;
% end
% 
% tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
% fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,tail_angle)
% 
% save('HALE_Validation_10ms.mat')

%% Perform full flight-dynamic simulation on trimmed/deformed aircraft
load('HALE_Validation_10ms_Rigid.mat')
COND.valSTAGGERSTEPS = 180;
COND.valGUSTAMP = 2;
COND.valGUSTL = 10;
COND.valGUSTSTART = 35;
FLAG.STIFFWING = 1;

TRIM.valCL = OUTP.vecCL(end);
TRIM.valCDI = OUTP.vecCDI(end);
TRIM.valE = OUTP.vecE(end);
TRIM.matVEHUVW = VEHI.matVEHUVW;

COND.valMAXTIME = 1500;
COND.valSTIFFSTEPS = 33;
COND.valSTARTFORCES = 1;

FLAG.FLIGHTDYN = 1;
FLAG.STEADY = 0;
FLAG.GUSTMODE = 0;

% De-rotate trim alpha from vehicle and adjust body frame velocities
% instead
% [ VEHI.matVEHUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, -COND.vecVEHALPHA, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
% [SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
%     INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
%     SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME);
% 
% [ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
%     SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, ~] ...
%     = fcnVLST2DVEPARAM(SURF.matNPDVE, SURF.matNPVLST);

% VEHI.matVEHUVW = [COND.vecVEHVINF*cosd(COND.vecVEHALPHA), 0, -COND.vecVEHVINF*sind(COND.vecVEHALPHA)];


[ ~, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, 90, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
[SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME] = fcnROTVEHICLEFLEX( SURF.matDVE, SURF.matNPDVE, SURF.matVLST, SURF.matCENTER,...
    INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH,...
    SURF.matNPVLST, INPU.vecVEHCG, SURF.matEALST, SURF.vecWINGCG, VEHI.vecPAYLCG, VEHI.vecFUSECG, VEHI.vecWINGCG(2:end,:), VEHI.vecBFRAME);

[ ~, ~, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
    SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, ~] ...
    = fcnVLST2DVEPARAM(SURF.matNPDVE, SURF.matNPVLST);

% VEHI.matVEHUVW(1) = VEHI.matVEHUVW(1)-0.1745;
VEHI.matGLOBUVW = VEHI.matVEHUVW;
VEHI.matVEHUVW = -VEHI.matVEHUVW;

TRIM.perturb(1:COND.valSTIFFSTEPS,4) = deg2rad(90);
VEHI.matROTMAT = [cos(pi-TRIM.perturb(1,4)) 0 -sin(pi-TRIM.perturb(1,4));...
                  0   1   0;...
                  sin(pi-TRIM.perturb(1,4)) 0 cos(pi-TRIM.perturb(1,4))]; % Rotation matrix from body frame to earth frame
VEHI.matGLOBUVW = (VEHI.matROTMAT*VEHI.matVEHUVW')'; % Vehicle velocity in the earth frame

[OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);    
