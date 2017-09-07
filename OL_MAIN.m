clc
clear
% delete('size.txt');

% warning off

% profile -memory on
disp('======================================================================================');
disp('  /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$$$$$$   /$$$$$$          /$$ /$$   /$$       ');
disp(' /$$__  $$| $$__  $$| $$_____/| $$__  $$ /$$__  $$        | $$|__/  | $$    ');
disp('| $$  \ $$| $$  \ $$| $$      | $$  \ $$| $$  \ $$        | $$ /$$ /$$$$$$    /$$$$$$ ');
disp('| $$  | $$| $$$$$$$/| $$$$$   | $$$$$$$/| $$$$$$$$ /$$$$$$| $$| $$|_  $$_/   /$$__  $$');
disp('| $$  | $$| $$____/ | $$__/   | $$__  $$| $$__  $$|______/| $$| $$  | $$    | $$$$$$$$');
disp('| $$  | $$| $$      | $$      | $$  \ $$| $$  | $$        | $$| $$  | $$ /$$| $$_____/');
disp('|  $$$$$$/| $$      | $$$$$$$$| $$  | $$| $$  | $$        | $$| $$  |  $$$$/|  $$$$$$$');
disp(' \______/ |__/      |________/|__/  |__/|__/  |__/        |__/|__/   \___/   \_______/');
disp('======================================================================================');
disp(' ');

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Reading in geometry
filename = 'inputs/simple_wing.vap';

[flagRELAX, flagSTEADY, matGEOM, valMAXTIME, valMINTIME, valDELTIME, valDELTAE, ...
    valDENSITY, valKINV, valVEHICLES, matVEHORIG, vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHROLL, ...
    vecVEHFPA, vecVEHTRK, ~, vecWINGTRI, vecWAKETRI, ~, vecAREA, vecSPAN, vecCMAC, ~, ...
    ~, vecSYM, vecN, vecM, ~, ~, ~, vecPANELWING, ...
    vecSURFACEVEHICLE, valPANELS, ~, vecROTORRPM, vecROTDIAM, matROTORHUB, matROTORAXIS, vecROTORBLADES, ~, vecPANELROTOR,...
    vecFTURB, vecFUSESECTIONS, matFGEOM, matSECTIONFUSELAGE, vecFUSEVEHICLE, matFUSEAXIS, matFUSEORIG, vecVEHRADIUS...
    ] = fcnXMLREAD(filename);

valMAXTIME = 0
flagRELAX = 0

vecM = [2 2]';
vecN = [10 10]';

vecWINGTRI(~isnan(vecWINGTRI)) = nan;
vecWAKETRI(~isnan(vecWAKETRI)) = nan;
flagTRI = 0;

flagPRINT   = 1;
flagPLOT    = 1;
flagCIRCPLOT = 1;
flagGIF = 0;
flagPREVIEW = 0;
flagPLOTWAKEVEL = 0;
flagPLOTUINF = 0;
flagVERBOSE = 0;

%% Discretizing geometry into DVEs
[matCENTER, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEROLL,...
    vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, matVLST, matNTVLST, matNPVLST, matDVE, valNELE,...
    matADJE, vecDVESYM, vecDVETIP, vecDVESURFACE, vecDVELE, vecDVETE, vecDVEPANEL, matPANELTE,...
    valWINGS,vecDVEVEHICLE, vecDVEWING, vecDVEROTOR, vecDVEROTORBLADE, matSURFACETYPE, vecROTORVEH, ...
    matFUSEGEOM, matVEHUVW, matVEHROT, matVEHROTRATE, matCIRORIG, vecVEHPITCH, vecVEHYAW,...
    matROTORHUBGLOB, matUINF, vecDVETRI, vecN, vecM, valWSIZE, valWSIZETRI] = fcnGEOM2DVE(matGEOM, ...
    matVEHORIG, vecWINGTRI, vecWAKETRI, vecN, vecM, vecPANELWING,...
    vecSYM, vecSURFACEVEHICLE, vecPANELROTOR, vecROTORBLADES, matROTORHUB, matROTORAXIS, matSECTIONFUSELAGE,...
    vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE, vecVEHVINF, vecVEHALPHA, vecVEHBETA, ...
    vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS, valVEHICLES, vecROTORRPM);

[hFig2] = fcnPLOTBODY(1, valNELE, matDVE, matVLST, matCENTER, []);

%% Add boundary conditions to D-Matrix

[matD] = fcnDWING_OL(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);
[matE] = fcnEWING(valNELE, matADJE, vecDVEHVCRD, vecDVELE, vecDVETE);
[matD] = fcnDEXPAND_OL(matD, matE, valNELE);

%% Add kinematic conditions to D-Matrix
[vecK] = fcnSINGFCT(valNELE, vecDVESURFACE, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON_OL(matD, valNELE, matDVE, matCENTER, matVLST, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecDVEHVCRD,vecSYM, vecDVELE);

%% Preparing to timestep
% Preallocating for a turbo-boost in performance
vecCL = nan(valMAXTIME,valVEHICLES);
vecCLF = nan(valMAXTIME,valVEHICLES);
vecCLI = nan(valMAXTIME,valVEHICLES);
vecCDI = nan(valMAXTIME,valVEHICLES);
vecE = nan(valMAXTIME,valVEHICLES);

% Initializing wake parameters
matWAKEGEOM = [];
matNPWAKEGEOM = [];
vecWDVEHVSPN = [];
vecWDVEHVCRD = [];
vecWDVEROLL = [];
vecWDVEPITCH = [];
vecWDVEYAW = [];
vecWDVELESWP = [];
vecWDVEMCSWP = [];
vecWDVETESWP = [];
vecWDVEAREA = [];
matWDVENORM = [];
matWVLST = [];
matWDVE = [];
valWNELE = 0;
matWCENTER = [];
matWCOEFF = [];
vecWK = [];
matWADJE = [];
vecWDVEPANEL = [];
valLENWADJE = 0;
vecWKGAM = [];
vecWDVESYM = [];
vecWDVETIP = [];
vecWDVESURFACE = [];
vecWDVETRI = [];
vecWDVEWING = [];
vecWKEGAM = [];

% Building wing resultant
[vecR] = fcnRWING_OL(valNELE, 0, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
    matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVETESWP, vecSYM, valWSIZE, vecDVELE, matVLST, matDVE);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED_OL(matD, vecR, valNELE);

for valTIMESTEP = 1:valMAXTIME
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
    
    [matUINF, matUINFTE, matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, ...
        matFUSEGEOM, matNEWWAKEPANEL] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, matVEHROTRATE, matCIRORIG, vecVEHRADIUS, ...
        valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE, matNPVLST, matFUSEGEOM, vecFUSEVEHICLE, ...
        matVEHROT, vecROTORVEH, matROTORHUBGLOB, matROTORHUB, matROTORAXIS, vecDVEROTOR, vecROTORRPM, matPANELTE);
    
    %% Generating new wake elements
    
    [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
        vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING, ...
        vecWDVELE, vecWDVETE, vecWKEGAM] = fcnCREATEWAKEROW_OL(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
        vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, ...
        vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVESURFACE, vecWDVEWING, flagSTEADY, valWSIZE, vecWKEGAM);
    
    %% Finding new wake coefficients
    
    [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
    [matWE] = fcnWEWAKE([1:valWNELE]', matWADJE, vecWDVEHVCRD, vecWDVELE, vecWDVETE, matCOEFF, vecDVEHVCRD, vecDVETE);
    [matWD] = fcnDEXPAND_OL(matWD, matWE, valWNELE);
    
    vecWR = [vecWR; zeros(size(matWE,1),1)];
    
    [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD_OL(matWD, vecWR, valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVEHVCRD(end-valWSIZE+1:end));
    
    %% Rebuilding and solving wing resultant
    [vecR] = fcnRWING_OL(valNELE, valTIMESTEP, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
        matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
        vecWDVETESWP, vecSYM, valWSIZE, vecDVELE, matVLST, matDVE);
    
    [matCOEFF] = fcnSOLVED_OL(matD, vecR, valNELE);
    
    %% Creating and solving WD-Matrix
    [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
    [matWE] = fcnWEWAKE([1:valWNELE]', matWADJE, vecWDVEHVCRD, vecWDVELE, vecWDVETE, matCOEFF, vecDVEHVCRD, vecDVETE);
    [matWD] = fcnDEXPAND_OL(matWD, matWE, valWNELE);
    
    vecWR = [vecWR; zeros(size(matWE,1),1)];
    
    [matWCOEFF] = fcnSOLVEWD_OL(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN, vecWDVEHVCRD);
    
    %% Relaxing wake
    if valTIMESTEP > 2 && flagRELAX == 1
        
        [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
            vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
            matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKE_OL(matUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
            matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVETE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
            vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
            vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING);
        
        % Creating and solving WD-Matrix
        [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
        [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
    end
    
    %% Forces
    %                 [vecCL(valTIMESTEP,:,ai), vecCLF(valTIMESTEP,:,ai), vecCLI(valTIMESTEP,:,ai), vecCDI(valTIMESTEP,:,ai), vecE(valTIMESTEP,:,ai), vecDVENFREE, vecDVENIND, ...
    %                     vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, ...
    %                     vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE, valWNELE, matWDVE, matWVLST, ...
    %                     matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, ...
    %                     vecSYM, vecDVETESWP, vecAREA, vecSPAN, [], vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEVEHICLE, valVEHICLES);
    
    %% Post-timestep outputs
    fcnPRINTOUT(flagPRINT, valTIMESTEP, valVEHICLES, vecCL, vecCDI)
    
    if flagGIF == 1 % Creating GIF (output to GIF/ folder by default)
        fcnGIF(flagVERBOSE, valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM, valWNELE, matWDVE, matWVLST, matWCENTER);
    end
    
end

fprintf('\n');

%% Plotting

if flagPLOT == 1
    fcnPLOTPKG(flagVERBOSE, flagPLOTWAKEVEL, 0, flagPLOTUINF, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM, valWNELE, matWDVE, matWVLST, matWCENTER, ...
        [], [], matUINF, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF);
    
    if flagCIRCPLOT == 1
        fcnPLOTCIRC_OL(valNELE, matDVE, matVLST, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF)
    end
end

