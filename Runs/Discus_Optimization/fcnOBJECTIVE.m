function [out] = fcnOBJECTIVE(design_var, N_bendstiff, N_torstiff, N_elasticaxis, N_massaxis, home_dir)

N_bendstiff = 19;
N_torstiff = 19;
N_elasticaxis = 19;
N_massaxis = 19;
home_dir = pwd;

cd(home_dir)
cd '../../'

% Write design variable history to file
fp2 = fopen('Optimization/dvhistory_sine.txt','at');
fprintf(fp2,'%g ', design_var);
fprintf(fp2,'\n');
fclose(fp2);


structure.matEIx(:,1) = design_var(1:N_bendstiff);
structure.matGJt(:,1) = design_var(N_bendstiff+1:N_bendstiff+N_torstiff);
structure.vecEA(:,1) = design_var(N_bendstiff+N_torstiff+1:N_bendstiff+N_torstiff+N_elasticaxis);
structure.vecCG(:,1) = design_var(N_bendstiff+N_torstiff+N_elasticaxis+1:N_bendstiff+N_torstiff+N_elasticaxis+N_massaxis);

baseline_file = 'inputs/Discus2c.vap';

trim_iter = 1;

VAP_IN = [];
TRIM = [];

% Initialize variables and read in geometry
[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(baseline_file,VAP_IN);

FLAG.OPT = 1; 

INPU.matEIx = structure.matEIx;
INPU.matGJt = structure.matGJt;
INPU.vecEA = structure.vecEA;
INPU.vecCG = structure.vecCG;

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

OUTP.TRIMFAIL = 0;
COND.valMAXTRIMITER = 50;

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

% tail_angle = rad2deg(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))./TRIM.tau;
fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)

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
        OUTP.aero_iter = 1;
        temp_def = [];
        temp_twist = [];
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

        if OUTP.aero_iter > COND.valMAXTRIMITER
            OUTP.TRIMFAIL = 1;
            break;
        end

        end

        if OUTP.TRIMFAIL == 0
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
        % Break out of function if stuck in trim loop
        if trim_iter > COND.valMAXTRIMITER || OUTP.TRIMFAIL == 1
            OUTP.TRIMFAIL = 1;
            break;
        end
    end
end

try
if OUTP.TRIMFAIL == 0
    OUTP.aero_iter = 0;
    fprintf('\nVehicle trimmed. AoA = %.2f deg., Elev. Angle = %.2f deg.\n\n',COND.vecVEHALPHA,SURF.vecELEVANGLE)

    %% Perform full flight-dynamic simulation on trimmed/deformed aircraft
    SURF.matBEAMACC = [];
    COND.valGUSTAMP = 1;
    COND.valGUSTSTART = 40;
    OUTP.valVS = COND.vecVEHVINF*sind(COND.vecVEHFPA);
    sink_rate = OUTP.valVS;

    SURF.matB = 0*[max(max(INPU.matEIx(:,1)))*8.333e-5; max(max(INPU.matGJt(:,1)))*1.6667e-4];

    COND.valMAXTIME = ceil((10*COND.valGUSTL + SURF.valTBOOM)/COND.vecVEHVINF/COND.valDELTIME + COND.valGUSTSTART);
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

    out = OUTP.vecZE_gain(end,1);
    
    [~,idxMAXdef] = max(OUTP.matDEFGLOB(:,end));
    [~,idxMAXtwist] = max(OUTP.matTWISTGLOB(:,end));
    [~,idxMINdef] = min(OUTP.matDEFGLOB(:,end));
    [~,idxMINtwist] = min(OUTP.matTWISTGLOB(:,end));
    struct = [OUTP.matDEFGLOB(idxMAXdef,:), OUTP.matDEFGLOB(idxMINdef,:), OUTP.matTWISTGLOB(idxMAXtwist,:), OUTP.matTWISTGLOB(idxMINtwist,:)];
    
else
    out = Inf;
    struct = Inf;
    sink_rate = Inf;
end
catch
    out = Inf;
    struct = Inf;
    sink_rate = Inf;
end
    

if nargin ~= 0
    fp2 = fopen('Optimization/opthistory_sine.txt','at');
    fprintf(fp2,'%g ', [out, sink_rate, design_var, struct]);
    fprintf(fp2,'\n');
    fclose(fp2);
end

cd(home_dir)


