function OUTP = fcnVAP_MAIN(filename, VAP_IN)

warning off
if nargin == 0
    VAP_MAIN;
    return
end

%% Reading in geometry
[FLAG, COND, VISC, INPU, VEHI] = fcnXMLREAD(filename, VAP_IN);

FLAG.PRINT = 1;
FLAG.PLOT = 1;
FLAG.VISCOUS = 1;
FLAG.CIRCPLOT = 0;
FLAG.GIF = 0;
FLAG.PREVIEW = 0;
FLAG.PLOTWAKEVEL = 0;
FLAG.PLOTUINF = 0;
FLAG.VERBOSE = 0;
FLAG.SAVETIMESTEP = 0;
FLAG.HOVERWAKE = 0;
FLAG.NACELLE = 0;
FLAG.GPU = 0;

% Initializing parameters to null/zero/nan
[WAKE, OUTP, INPU, SURF] = fcnINITIALIZE(COND, INPU);

if FLAG.PRINT == 1
    disp('============================================================================');
    disp('                  /$$    /$$  /$$$$$$  /$$$$$$$         /$$$$$$     /$$$$$$$') 
    disp('+---------------+| $$   | $$ /$$__  $$| $$__  $$       /$$__  $$   | $$____/') ;
    disp('| RYERSON       || $$   | $$| $$  \ $$| $$  \ $$      |__/  \ $$   | $$      ');
    disp('| APPLIED       ||  $$ / $$/| $$$$$$$$| $$$$$$$/         /$$$$$/   | $$$$$$$ ');
    disp('| AERODYNAMICS  | \  $$ $$/ | $$__  $$| $$____/         |___  $$   |_____  $$');
    disp('| LABORATORY OF |  \  $$$/  | $$  | $$| $$             /$$  \ $$    /$$  \ $$');
    disp('| FLIGHT        |   \  $/   | $$  | $$| $$            |  $$$$$$//$$|  $$$$$$/');
    disp('+---------------+    \_/    |__/  |__/|__/             \______/|__/ \______/ ');
    disp('============================================================================');
    disp(' ');                                                   
end

% Setting up timestep saving feature
if FLAG.SAVETIMESTEP == 1
    if exist('timesteps/') ~= 7; mkdir('timesteps'); end
    if isfield(VAP_IN,'TimestepName')
        timestep_folder = strcat('timesteps/',VAP_IN.TimestepName,'/');
    else
        timestep_folder = ['timesteps/',regexprep(filename,{'inputs/', '.vap'}, ''), '_(', datestr(now, 'dd_mm_yyyy HH_MM_SS_FFF'),')/'];
    end
    mkdir(timestep_folder); 
end

% Check if the files required by the viscous calculations exist
[FLAG] = fcnVISCOUSFILECHECK(FLAG, VISC);

%% Discretizing geometry into DVEs
% Adding collective pitch to the propeller/rotor
if ~isempty(COND.vecCOLLECTIVE)
    INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) = INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) + repmat(reshape(COND.vecCOLLECTIVE(INPU.vecPANELROTOR(INPU.vecPANELROTOR > 0), 1),1,1,[]),2,1,1);
end
[INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, OUTP, SURF);
% fcnPLOTPKG([], FLAG, SURF, VISC, WAKE, COND, INPU)
%% Advance Ratio
MISC.vecROTORJ = [];
for jj = 1:length(COND.vecROTORRPM)
    MISC.vecROTORJ(jj) = (COND.vecVEHVINF(VEHI.vecROTORVEH(jj))*60)./(abs(COND.vecROTORRPM(jj)).*INPU.vecROTDIAM(jj));
end

%% Add boundary conditions to D-Matrix
[matD] = fcnDWING(SURF, INPU);

%% Add kinematic conditions to D-Matrix
[SURF.vecK] = fcnSINGFCT(SURF.valNELE, SURF.vecDVESURFACE, SURF.vecDVETIP, SURF.vecDVEHVSPN);
[matD] = fcnKINCON(matD, SURF, INPU, FLAG);

%% Preparing to timestep
% Building wing resultant
[vecR] = fcnRWING(0, SURF, WAKE, FLAG);

% Solving for wing coefficients
[SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
SURF.matNPDVE = SURF.matDVE;

% Computing structure distributions if data exists
try 
    [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF); 
    FLAG.STRUCTURE = 1; % Create flag if structure data exists
catch
    FLAG.STRUCTURE = 0; 
end

n = 1;
valGUSTTIME = 1;
SURF.gust_vel_old = zeros(SURF.valNELE,1);

%% Timestepping
for valTIMESTEP = 1:COND.valMAXTIME
    %% Timestep to solution
    %   Move wing
    %   Generate new wake elements
    %   Create and solve WD-Matrix for new elements
    %   Solve wing D-Matrix with wake-induced velocities
    %   Solve entire WD-Matrix
    %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
    %   Calculate surface normal forces
    %   Calculate DVE normal forces
    %   Calculate induced drag
    %   Calculate cn, cl, cy, cdi
    %   Calculate viscous effects
    
    %% Moving the vehicles
    
    % Bend wing if applicable, else move wing normally
    if FLAG.STRUCTURE == 1
        [SURF, INPU, COND, MISC, VISC, OUTP, valGUSTTIME, n] = fcnMOVESTRUCTURE(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valGUSTTIME, valTIMESTEP, n);
    else
        [SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);
    end
    
    %% Generating new wake elements
    [INPU, COND, MISC, VISC, WAKE, VEHI, SURF] = fcnCREATEWAKEROW(FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF);
    
    if FLAG.PREVIEW ~= 1
        %% Creating and solving WD-Matrix for latest row of wake elements
        % We need to grab from WAKE.matWADJE only the values we need for this latest row of wake DVEs
        idx = sparse(sum(ismember(WAKE.matWADJE,[((WAKE.valWNELE - WAKE.valWSIZE) + 1):WAKE.valWNELE]'),2)>0 & (WAKE.matWADJE(:,2) == 4 | WAKE.matWADJE(:,2) == 2));
        temp_WADJE = [WAKE.matWADJE(idx,1) - (valTIMESTEP-1)*WAKE.valWSIZE WAKE.matWADJE(idx,2) WAKE.matWADJE(idx,3) - (valTIMESTEP-1)*WAKE.valWSIZE];
        
        [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWSIZE]', temp_WADJE, WAKE.vecWDVEHVSPN(end-WAKE.valWSIZE+1:end), WAKE.vecWDVESYM(end-WAKE.valWSIZE+1:end), WAKE.vecWDVETIP(end-WAKE.valWSIZE+1:end), WAKE.vecWKGAM(end-WAKE.valWSIZE+1:end), INPU.vecN);
        [WAKE.matWCOEFF(end-WAKE.valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWSIZE, WAKE.vecWKGAM(end-WAKE.valWSIZE+1:end), WAKE.vecWDVEHVSPN(end-WAKE.valWSIZE+1:end));
        
        %% Rebuilding and solving wing resultant
        [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
        [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
        
        %% Creating and solving WD-Matrix
        [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
        [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
        
        %% Relaxing wake
        if valTIMESTEP > 2 && FLAG.RELAX == 1
            old_span = WAKE.vecWDVEHVSPN;
            WAKE = fcnRELAXWAKE(valTIMESTEP, SURF, WAKE, COND, FLAG, INPU);
            WAKE.matWCOEFF(:,2:3) = WAKE.matWCOEFF(:,2:3).*[old_span./WAKE.vecWDVEHVSPN (old_span./WAKE.vecWDVEHVSPN).^2];
        end
        
        %% Forces
        if valTIMESTEP >= COND.valSTARTFORCES
            [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP);
        end
        
        if FLAG.SAVETIMESTEP == 1
            save([timestep_folder, 'timestep_', num2str(valTIMESTEP), '.mat'], 'filename','valTIMESTEP','INPU','COND','MISC','WAKE','VEHI','SURF','OUTP');
        end
    end
    
    %% Post-timestep outputs
    if FLAG.PRINT == 1
        fcnPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.vecCT, MISC.vecROTORJ, VEHI.vecROTORVEH, 1)
    end
    
    if FLAG.GIF == 1 % Creating GIF (output to GIF/ folder by default)
       fcnGIF(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU, 1)
    end
end

[OUTP] = fcnOUTPUT(COND, FLAG, INPU, SURF, OUTP, valTIMESTEP);

if FLAG.PRINT == 1 && FLAG.PREVIEW == 0
    fprintf('VISCOUS CORRECTIONS => CLv = %0.4f \tCD = %0.4f \n', OUTP.vecCLv(end,:), OUTP.vecCD(end,:))
    fprintf('\n');
end

%% Plotting
if FLAG.PLOT == 1
    fcnPLOTPKG(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU)
end