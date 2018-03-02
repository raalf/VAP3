clc
clear
warning off

try
profile -memory on
profile on
end

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

%% Best Practices
% 1. Define wing from one wingtip to another in one direction

%% Reading in geometry
% filename = 'inputs/J_COLE_BASELINE.vap';
% filename = 'inputs/J_COLE_X57_CRUISE_PROP.vap'
filename = 'inputs/QuadRotor.vap'
% filename = 'inputs/simple-wing.vap'

[FLAG, COND, VISC, INPU, VEHI] = fcnXMLREAD(filename);
FLAG.RELAX = 1
COND.valMAXTIME = 2

COND.vecWINGTRI(~isnan(COND.vecWINGTRI)) = nan;
COND.vecWAKETRI(~isnan(COND.vecWAKETRI)) = nan;
FLAG.TRI = 0;
FLAG.GPU = 0;

FLAG.PRINT   = 1;
FLAG.PLOT    = 1;
FLAG.CIRCPLOT = 0;
FLAG.GIF = 0;
FLAG.PREVIEW = 0;
FLAG.PLOTWAKEVEL = 0;
FLAG.PLOTUINF = 0;
FLAG.VERBOSE = 0;
FLAG.VISCOUS = 0;

% Initializing parameters to null/zero/nan
[WAKE, OUTP, INPU] = fcnINITIALIZE(COND, INPU);

%% Discretizing geometry into DVEs
% Adding collective pitch to the propeller/rotor
%     tINPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) = INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) + repmat(reshape(vecCOLLECTIVE(INPU.vecPANELROTOR(INPU.vecPANELROTOR > 0), 1),1,1,[]),2,1,1);
[INPU, COND, MISC, VISC, WAKE, VEHI, SURF] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE);

%% Advance Ratio
MISC.vecROTORJ = [];
for jj = 1:length(COND.vecROTORRPM)
    MISC.vecROTORJ(jj) = (COND.vecVEHVINF(VEHI.vecROTORVEH(jj))*60)./(abs(COND.vecROTORRPM(jj)).*INPU.vecROTDIAM(jj));
end

%% Add boundary conditions to D-Matrix
[matD] = fcnDWING(SURF.valNELE, SURF.matADJE, SURF.vecDVEHVSPN, SURF.vecDVESYM, SURF.vecDVETIP, INPU.vecN);

%% Add kinematic conditions to D-Matrix
[SURF.vecK] = fcnSINGFCT(SURF.valNELE, SURF.vecDVESURFACE, SURF.vecDVETIP, SURF.vecDVEHVSPN);
[matD] = fcnKINCON(matD, SURF.valNELE, SURF.matDVE, SURF.matCENTER, SURF.matVLST, SURF.matDVENORM, SURF.vecK, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, SURF.vecDVELESWP, SURF.vecDVETESWP, SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, INPU.vecSYM, FLAG.GPU);

%% Preparing to timestep
% Building wing resultant
[vecR] = fcnRWING(SURF.valNELE, 0, SURF.matCENTER, SURF.matDVENORM, SURF.matUINF, WAKE.valWNELE, WAKE.matWDVE, ...
    WAKE.matWVLST, WAKE.matWCOEFF, WAKE.vecWK, WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD,WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, ...
    WAKE.vecWDVETESWP, SURF.vecDVESYM, WAKE.valWSIZE, FLAG.TRI, FLAG.STEADY);

% Solving for wing coefficients
[SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);

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
    
    [SURF.matUINF, SURF.matUINFTE, INPU.matVEHORIG, SURF.matVLST, SURF.matCENTER, MISC.matNEWWAKE, MISC.matNPNEWWAKE, ...
        VISC.matFUSEGEOM, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, SURF.matDVENORM, SURF.matNTVLST, SURF.matUINFROT] = fcnMOVESURFACE(INPU.matVEHORIG, VEHI.matVEHUVW, VEHI.matVEHROTRATE, MISC.matCIRORIG, VEHI.vecVEHRADIUS, ...
        COND.valDELTIME, SURF.matVLST, SURF.matCENTER, SURF.matDVE, SURF.vecDVEVEHICLE, SURF.vecDVETE, VISC.matFUSEGEOM, VISC.vecFUSEVEHICLE, ...
        VEHI.matVEHROT, VEHI.vecROTORVEH, INPU.matROTORHUB, INPU.matROTORAXIS, SURF.vecDVEROTOR, COND.vecROTORRPM, SURF.matPANELTE, SURF.matNTVLST);
    
    %% Generating new wake elements
    
    idxtri = SURF.vecDVETRI == 1;
    widxtri = WAKE.vecWDVETRI == 1;
    
    [INPU, COND, MISC, VISC, WAKE, VEHI, SURF] = fcnCREATEWAKEROW(FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF);
    
    if FLAG.PREVIEW ~= 1
        %% Creating and solving WD-Matrix for latest row of wake elements
        % We need to grab from WAKE.matWADJE only the values we need for this latest row of wake DVEs
        idx = sparse(sum(ismember(WAKE.matWADJE,[((WAKE.valWNELE - WAKE.valWSIZE) + 1):WAKE.valWNELE]'),2)>0 & (WAKE.matWADJE(:,2) == 4 | WAKE.matWADJE(:,2) == 2));
        temp_WADJE = [WAKE.matWADJE(idx,1) - (valTIMESTEP-1)*WAKE.valWSIZE WAKE.matWADJE(idx,2) WAKE.matWADJE(idx,3) - (valTIMESTEP-1)*WAKE.valWSIZE];
        
        [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWSIZE]', temp_WADJE, WAKE.vecWDVEHVSPN(end-WAKE.valWSIZE+1:end), WAKE.vecWDVESYM(end-WAKE.valWSIZE+1:end), WAKE.vecWDVETIP(end-WAKE.valWSIZE+1:end), WAKE.vecWKGAM(end-WAKE.valWSIZE+1:end), INPU.vecN);
        [WAKE.matWCOEFF(end-WAKE.valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWSIZE, WAKE.vecWKGAM(end-WAKE.valWSIZE+1:end), WAKE.vecWDVEHVSPN(end-WAKE.valWSIZE+1:end));
        
        %% Rebuilding and solving wing resultant
        [vecR] = fcnRWING(SURF.valNELE, valTIMESTEP, SURF.matCENTER, SURF.matDVENORM, SURF.matUINF, WAKE.valWNELE, WAKE.matWDVE, ...
            WAKE.matWVLST, WAKE.matWCOEFF, WAKE.vecWK, WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD,WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, ...
            WAKE.vecWDVETESWP, SURF.vecDVESYM, WAKE.valWSIZE, FLAG.TRI, FLAG.STEADY, FLAG.GPU);
        
        [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
        
        %% Creating and solving WD-Matrix
        [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
        [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
        
        %% Relaxing wake
        if valTIMESTEP > 2 && FLAG.RELAX == 1
            WAKE = fcnRELAXWAKE(valTIMESTEP, SURF, WAKE, COND, FLAG);

            % Creating and solving WD-Matrix
            [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
            [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
        end
        
        %% Forces
        [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnFORCES(valTIMESTEP, FLAG, INPU, COND, MISC, VISC, WAKE, VEHI, SURF);
  
    end
    
    %% Post-timestep outputs
    if FLAG.PRINT == 1
        fcnPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.vecCTCONV, MISC.vecROTORJ, VEHI.vecROTORVEH, 1)
    end
    
    if FLAG.GIF == 1 % Creating GIF (output to GIF/ folder by default)
        fcnGIF(FLAG.VERBOSE, valTIMESTEP, SURF.valNELE, SURF.matDVE, SURF.matVLST, SURF.matCENTER, VISC.matFUSEGEOM, WAKE.valWNELE, WAKE.matWDVE, WAKE.matWVLST, WAKE.matWCENTER, WAKE.vecWPLOTSURF, 1);
    end
    
    if FLAG.PREVIEW ~= 1 && max(SURF.vecDVEROTOR) > 0
        temp_cdi = fcnTIMEAVERAGE(OUTP.vecCDI(:,end), COND.vecROTORRPM, COND.valDELTIME);
    end
end

if FLAG.PRINT == 1 && FLAG.PREVIEW == 0
    fprintf('VISCOUS CORRECTIONS => CLv = %0.4f \tCD = %0.4f \n', OUTP.vecCLv(end,:), OUTP.vecCD(end,:))
    fprintf('\n');
end

%% Plotting

fcnPLOTPKG(FLAG, SURF, VISC, WAKE, COND)

profile report
profile off
