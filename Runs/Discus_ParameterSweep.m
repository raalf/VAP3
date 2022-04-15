clc
clear
warning off

cores = 64;
parpool(cores,'IdleTimeout',800)

cd '..'

addpath('Flight Dynamics')

filename = 'inputs/Discus2c_Rect.vap';

delete Optimization/paramhistory.txt

delete Optimization/dvparamhistory.txt

matEIx = [100000 250000 500000 1000000 1500000];
matGJt = [100000 250000 500000 1000000 1500000];
EA = [0.2 0.3 0.4 0.5 0.6];
CG = [0.2 0.3 0.4 0.5 0.6];

param_sweep = combvec(matEIx, matGJt, EA, CG);

parfor kk = 1:size(param_sweep,2)
% for kk = 1:size(param_sweep,2)
    
fp3 = fopen('Optimization/dvparamhistory.txt','at');
fprintf(fp3,'%g ', param_sweep(:,kk)');
fprintf(fp3,'\n');
fclose(fp3);

temp_def = [];
temp_twist = [];
trim_def = [];
trim_twist = [];
trim_alpha = [];
temp = [];
trim_iter = 1;

VAP_IN = [];
TRIM = [];
OUTP.TRIMFAIL = 0;

% Initialize variables and read in geometry
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN);

FLAG.OPT = 2;
COND.valMAXTRIMITER = 50;

INPU.matEIx_param = param_sweep(1,kk);
INPU.matGJt_param = param_sweep(2,kk);
INPU.vecEA_param = param_sweep(3,kk);
INPU.vecCG_param = param_sweep(4,kk);

[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP);

CGX_loc = 0.40;
Xnew = (CGX_loc*COND.vecVEHWEIGHT/9.81 - VEHI.vecWINGMASS(1)*VEHI.vecWINGCG(1,1))/sum(VEHI.vecFUSEMASS,1);
[INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
dX = INPU.vecVEHCG(1) - CGX_loc;

VEHI.vecFUSELOC = [Xnew - VEHI.vecFUSEL/2,0,0];
VEHI.vecFUSEMASS = 280.4158;
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
FLAG.VISCOUS = 1;
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
if OUTP.TRIMFAIL == 0
while max(abs(tol)) > 0.01
    COND.valSTIFFSTEPS = COND.valMAXTIME - 1;
    
    COND.valSTARTFORCES = COND.valMAXTIME - 2; % Compute forces for last two timesteps
    
    tol_aero = 100;
    tol_aero2 = 100;
    tol_aero1 = 100;
    SURF.matB = [max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];
    
    %% Aeroelastic convergence loop
    % Compute static aeroelastic configuration based on trimmed aircraft loads
    while max(abs(tol_aero)) > 1e-3 && OUTP.TRIMFAIL == 0
        
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
    
    % Break out of function if stuck in trim loop
    if trim_iter > COND.valMAXTRIMITER || OUTP.TRIMFAIL == 1
        OUTP.TRIMFAIL = 1;
        break;
    end
        
end

OUTP.aero_iter = 0;
tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)
end
end

try
if OUTP.TRIMFAIL == 0
    temp.OUTP = OUTP;
    temp.COND = COND;
    temp.INPU = INPU;
    temp.FLAG = FLAG;
    temp.MISC = MISC;
    temp.SURF = SURF;
    temp.TRIM = TRIM;
    temp.VEHI = VEHI;
    temp.VISC = VISC;
    temp.WAKE = WAKE;
    
    %% Perform full flight-dynamic simulation on trimmed/deformed aircraft
    SURF.matBEAMACC = [];
    COND.valGUSTAMP = 1;
    COND.valGUSTSTART = 40;
    OUTP.valVS = COND.vecVEHVINF*sind(COND.vecVEHFPA);
    sink_rate = OUTP.valVS;

    SURF.matB = 0*[max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];

    COND.valMAXTIME = ceil((COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART);
    COND.valSTIFFSTEPS = 15;
    COND.valSTARTFORCES = 1;
    FLAG.FLIGHTDYN = 1;
    FLAG.STATICAERO = 0;
    FLAG.STEADY = 0;
    FLAG.RELAX = 0;
    FLAG.GUSTMODE = 1;
    FLAG.SAVETIMESTEP = 0;

    VEHI.vecVEHDYN(1:COND.valSTIFFSTEPS,4) = deg2rad(COND.vecVEHPITCH);

    [VEHI.matVEHUVW] = fcnGLOBSTAR(VEHI.matGLOBUVW, 0, pi+deg2rad(COND.vecVEHPITCH), 0);

    COND.start_loc = repmat([-COND.valGUSTSTART*COND.valDELTIME*COND.vecVEHVINF,0,0],size(SURF.matCENTER,1),1, size(SURF.matCENTER,3)); % Location (in meters) in global frame where gust starts

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, 1);

    gain(kk,1) = OUTP.vecZE_gain(end);

else
    gain(kk,1) = Inf;
    sink_rate = Inf;
end
catch
    gain(kk,1) = Inf;
    sink_rate = Inf;
end

fp2 = fopen('Optimization/paramhistory_rect_sine.txt','at');
fprintf(fp2,'%g ', [gain(kk,1), sink_rate, param_sweep(:,kk)']);
fprintf(fp2,'\n');
fclose(fp2);

end
