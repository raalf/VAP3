function OUTP = fcnVAP_MAIN(filename, VAP_IN)

warning off
if nargin == 0
    VAP_MAIN;
    return
end

%% Reading in geometry
[FLAG, COND, VISC, INPU, VEHI] = fcnXMLREAD(filename, VAP_IN);

COND.vecWINGTRI(~isnan(COND.vecWINGTRI)) = nan;
COND.vecWAKETRI(~isnan(COND.vecWAKETRI)) = nan;
FLAG.TRI = 0;
FLAG.GPU = 0;

FLAG.PRINT = 0;
FLAG.PLOT = 0;
FLAG.VISCOUS = 1;
FLAG.CIRCPLOT = 0;
FLAG.GIF = 0;
FLAG.PREVIEW = 0;
FLAG.PLOTWAKEVEL = 0;
FLAG.PLOTUINF = 0;
FLAG.VERBOSE = 0;
FLAG.SAVETIMESTEP = 0;
FLAG.NACELLE = 1;

% Initializing parameters to null/zero/nan
[WAKE, OUTP, INPU, SURF] = fcnINITIALIZE(COND, INPU);

if FLAG.PRINT == 1
    disp('============================================================================');
    disp('                  /$$    /$$  /$$$$$$  /$$$$$$$         /$$$$$$      /$$');
    disp('+---------------+| $$   | $$ /$$__  $$| $$__  $$       /$$__  $$    /$$$$');
    disp('| RYERSON       || $$   | $$| $$  \ $$| $$  \ $$      |__/  \ $$   |_  $$');
    disp('| APPLIED       ||  $$ / $$/| $$$$$$$$| $$$$$$$/         /$$$$$/     | $$');
    disp('| AERODYNAMICS  | \  $$ $$/ | $$__  $$| $$____/         |___  $$     | $$');
    disp('| LABORATORY OF |  \  $$$/  | $$  | $$| $$             /$$  \ $$     | $$');
    disp('| FLIGHT        |   \  $/   | $$  | $$| $$            |  $$$$$$//$$ /$$$$$$');
    disp('+---------------+    \_/    |__/  |__/|__/             \______/|__/|______/');
    disp('============================================================================');
    disp(' ');
end

% Setting up timestep saving feature
if FLAG.SAVETIMESTEP == 1
    if exist('timesteps/') ~= 7; mkdir(timesteps); end
    timestep_folder = ['timesteps/',regexprep(filename,{'inputs/', '.vap'}, ''), '_(', datestr(now, 'dd_mm_yyyy HH_MM_SS_FFF'),')/'];
    mkdir(timestep_folder);
end

% Check if the files required by the viscous calculations exist
[FLAG] = fcnVISCOUSFILECHECK(FLAG, VISC);

%% Discretizing geometry into DVEs
% Adding collective pitch to the propeller/rotor
if ~isempty(COND.vecCOLLECTIVE)
    INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) = INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) + repmat(reshape(COND.vecCOLLECTIVE(INPU.vecPANELROTOR(INPU.vecPANELROTOR > 0), 1),1,1,[]),2,1,1);
end
[INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, FLAG, OUTP);
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
        if valTIMESTEP <= COND.valSTIFFSTEPS || FLAG.STIFFWING == 1
            
            [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, valTIMESTEP);
            
            if FLAG.STIFFWING == 2
                [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF);
            end
            
        elseif valTIMESTEP == n*COND.valSTIFFSTEPS + 1 || valGUSTTIME > 1
            [COND, INPU, OUTP, MISC, SURF, valGUSTTIME] = fcnFLEXWING(INPU, COND, SURF, OUTP, FLAG, MISC, valGUSTTIME, valTIMESTEP);
            n = n + 1;
        else
            [SURF, INPU, MISC, VISC, OUTP] = fcnSTIFFWING_STATIC(INPU, VEHI, MISC, COND, SURF, VISC, OUTP, valTIMESTEP);
        end
    else
        [SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);
    end
    
    % Update structure location after moving wing
    try [SURF] = fcnWINGSTRUCTGEOM(SURF, INPU); 
    catch 
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
            WAKE = fcnRELAXWAKE(valTIMESTEP, SURF, WAKE, COND, FLAG, INPU);
            
            % Creating and solving WD-Matrix
            [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
            [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
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
        fcnGIF(FLAG.VERBOSE, valTIMESTEP, SURF.valNELE, SURF.matDVE, SURF.matVLST, SURF.matCENTER, VISC.matFUSEGEOM, WAKE.valWNELE, WAKE.matWDVE, WAKE.matWVLST, WAKE.matWCENTER, WAKE.vecWPLOTSURF, 1);
    end
    
end

if FLAG.PREVIEW ~= 1 && max(SURF.vecDVEROTOR) > 0 && ~isempty(valTIMESTEP)
    % Time averaged lift
    OUTP.vecCL_AVG = fcnTIMEAVERAGE(OUTP.vecCLv, COND.vecROTORRPM, COND.valDELTIME);
    
    % Time averaged drags (total, induced, profile)
    OUTP.vecCD_AVG = fcnTIMEAVERAGE(OUTP.vecCD, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCDI_AVG = fcnTIMEAVERAGE(OUTP.vecCDI, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCDP_AVG = fcnTIMEAVERAGE(OUTP.vecCD - OUTP.vecCDI, COND.vecROTORRPM, COND.valDELTIME);
    
%     for i = 1:max(SURF.vecDVEWING)
%        OUTP.WING(i).vecLDIST(~any(OUTP.WING(i).vecLDIST, 2), :) = [];
%        OUTP.WING(i).vecLDIST_AVG = fcnTIMEAVERAGE(OUTP.WING(i).vecLDIST, COND.vecROTORRPM, COND.valDELTIME);
%        OUTP.WING(i).vecDPDIST(~any(OUTP.WING(i).vecDPDIST, 2), :) = [];
%        OUTP.WING(i).vecDPDIST_AVG = fcnTIMEAVERAGE(OUTP.WING(i).vecDPDIST, COND.vecROTORRPM, COND.valDELTIME);
%        OUTP.WING(i).vecDIDIST(~any(OUTP.WING(i).vecDIDIST, 2), :) = [];
%        OUTP.WING(i).vecDIDIST_AVG = fcnTIMEAVERAGE(OUTP.WING(i).vecDIDIST, COND.vecROTORRPM, COND.valDELTIME);
%     end
    
    OUTP.vecCT_AVG = fcnTIMEAVERAGE(OUTP.vecCT, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCFx_AVG = fcnTIMEAVERAGE(OUTP.vecCFx, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCFy_AVG = fcnTIMEAVERAGE(OUTP.vecCFy, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCMx_AVG = fcnTIMEAVERAGE(OUTP.vecCMx, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCMy_AVG = fcnTIMEAVERAGE(OUTP.vecCMy, COND.vecROTORRPM, COND.valDELTIME);
    
    % Time averaged propeller powers
    OUTP.vecCP_AVG = fcnTIMEAVERAGE(OUTP.vecCP, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCPI_AVG = fcnTIMEAVERAGE(OUTP.vecCPI, COND.vecROTORRPM, COND.valDELTIME);
    OUTP.vecCPP_AVG = fcnTIMEAVERAGE(OUTP.vecCP - OUTP.vecCPI, COND.vecROTORRPM, COND.valDELTIME);

    for i = 1:max(SURF.vecDVEROTOR)
        OUTP.ROTOR(i).vecTHRUSTDIST(any(~any(OUTP.ROTOR(i).vecTHRUSTDIST, 2), 3), :) = [];
        OUTP.ROTOR(i).vecTHRUSTDIST_AVG = fcnTIMEAVERAGE(OUTP.ROTOR(i).vecTHRUSTDIST, COND.vecROTORRPM, COND.valDELTIME);
        OUTP.ROTOR(i).vecTORQUEDIST(any(~any(OUTP.ROTOR(i).vecTORQUEDIST, 2), 3), :) = [];
        OUTP.ROTOR(i).vecTORQUEDIST_AVG = fcnTIMEAVERAGE(OUTP.ROTOR(i).vecTORQUEDIST, COND.vecROTORRPM, COND.valDELTIME);
    end
    
elseif FLAG.PREVIEW ~= 1 && ~any(SURF.vecDVEROTOR)
    OUTP.vecCLv_AVG = OUTP.vecCLv(end);
    OUTP.vecCD_AVG = OUTP.vecCD(end);
    OUTP.vecCDI_AVG = OUTP.vecCDI(end);
    OUTP.vecCDP_AVG = OUTP.vecCD(end) - OUTP.vecCDI(end);
end

if FLAG.PRINT == 1 && FLAG.PREVIEW == 0
    fprintf('VISCOUS CORRECTIONS => CLv = %0.4f \tCD = %0.4f \n', OUTP.vecCLv(end,:), OUTP.vecCD(end,:))
    fprintf('\n');
end

%% Plotting
if FLAG.PLOT == 1
    fcnPLOTPKG(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU)
end

OUTP.vecVEHALPHA = COND.vecVEHALPHA;
OUTP.vecCOLLECTIVE = COND.vecCOLLECTIVE;
OUTP.vecROTDIAM = INPU.vecROTDIAM;
OUTP.vecVEHWEIGHT = COND.vecVEHWEIGHT;
OUTP.vecROTORRPM = COND.vecROTORRPM;
OUTP.vecDVEAREA = SURF.vecDVEAREA;
OUTP.valAREA = INPU.vecAREA;
OUTP.matGEOM = INPU.matGEOM;
OUTP.valDENSITY = COND.valDENSITY;