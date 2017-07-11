% clc
clear
% delete('size.txt');

warning off

% profile -memory on

disp('=============================================================================');
disp('                  /$$    /$$  /$$$$$$  /$$$$$$$         /$$$$$$      /$$$$$$ ');
disp('+---------------+| $$   | $$ /$$__  $$| $$__  $$       /$$__  $$    /$$$_  $$');
disp('| RYERSON       || $$   | $$| $$  \ $$| $$  \ $$      |__/  \ $$   | $$$$\ $$');
disp('| APPLIED       ||  $$ / $$/| $$$$$$$$| $$$$$$$/         /$$$$$/   | $$ $$ $$');
disp('| AERODYNAMICS  | \  $$ $$/ | $$__  $$| $$____/         |___  $$   | $$\ $$$$');
disp('| LABORATORY OF |  \  $$$/  | $$  | $$| $$             /$$  \ $$   | $$ \ $$$');
disp('| FLIGHT        |   \  $/   | $$  | $$| $$            |  $$$$$$//$$|  $$$$$$/');
disp('+---------------+    \_/    |__/  |__/|__/             \______/|__/ \______/');
disp('=============================================================================');
disp(' ');

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Reading in geometry

% filename = 'inputs/2MotorGliders.vap';
% filename = 'inputs/QuadRotor.vap';
% filename = 'inputs/TMotor.vap';
% filename = 'inputs/StandardCirrusSym.vap';
% filename = 'inputs/StandardCirrus.vap';
% filename = 'inputs/StandardCirrusTail2.vap';
% filename = 'inputs/XMLtest.vap';
filename = 'inputs/twoVehicles.vap';

[flagRELAX, flagSTEADY, matGEOM, valMAXTIME, valMINTIME, valDELTIME, valDELTAE, ...
    valDENSITY, valKINV, valVEHICLES, matVEHORIG, vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHROLL, ...
    vecVEHFPA, vecVEHTRK, ~, vecWINGTRI, vecWAKETRI, ~, vecAREA, vecSPAN, vecCMAC, ~, ...
    ~, vecSYM, vecN, vecM, ~, ~, ~, ~, ...
    vecSURFACEVEHICLE, valPANELS, ~, vecROTORRPM, vecROTDIAM, matROTORHUB, matROTORAXIS, vecROTORBLADES, ~, vecROTOR,...
    vecFTURB, vecFUSESECTIONS, matFGEOM, matSECTIONFUSELAGE, vecFUSEVEHICLE, matFUSEAXIS, matFUSEORIG, vecVEHRADIUS...
    ] = fcnXMLREAD(filename);

% For debugging:
valMAXTIME = 100
% vecVEHFPA = 0
% vecVEHTRK = 0


flagRELAX = 0;

flagTRI = 0;

flagPRINT   = 1;
flagPLOT    = 1;
flagCIRCPLOT = 0;
flagGIF = 0;
flagPREVIEW = 1;
flagPLOTWAKEVEL = 0;
flagPLOTUINF = 0;
flagVERBOSE = 0;


%%

% tranlsate matGEOM to vehicle origin
matGEOM(:,1:3,:) = matGEOM(:,1:3,:)+permute(reshape(matVEHORIG(matGEOM(:,6,:),:)',3,2,[]),[2,1,3]);

[matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVESURFACE, vecDVELE, vecDVETE, vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);


% % Identifying which DVEs belong to which vehicle, as well as which type of lifting surface they belong to (wing or rotor)
vecDVEVEHICLE = vecSURFACEVEHICLE(vecDVESURFACE);
vecDVEWING = vecDVESURFACE;

vecDVEROTOR = vecROTOR(vecDVEPANEL); % Alton-Y
vecDVEROTORBLADE = vecDVEROTOR; % Current rotor DVEs are for Blade 1 (they are duplicated to Blade 2, 3, etc etc below)
idx_rotor = vecDVEROTOR>0; % Alton-Y
vecDVEWING(idx_rotor) = 0;

matSURFACETYPE = zeros(size(unique(vecDVESURFACE),1),2);
matSURFACETYPE(nonzeros(unique(vecDVEWING)),1) = nonzeros(unique(vecDVEWING));
matSURFACETYPE(nonzeros(unique(vecDVESURFACE(idx_rotor))),2) = nonzeros(unique(vecDVEROTOR));


% Identifying which ROTOR belongs to which vehicle.
vecROTORVEH = vecSURFACEVEHICLE(matSURFACETYPE(:,2)~=0);

% Duplicate Blades in a Rotor
[ matVLST0, matNPVLST0, matCENTER0, matDVE, matADJE, vecDVEVEHICLE, ...
    vecDVEWING, vecDVEROTOR, matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, ...
    vecDVETIP, vecDVELE, vecDVETE, vecDVEROTORBLADE, vecDVESYM, ...
    valNELE ] = fcnDUPBLADE( vecROTORVEH, vecDVEROTOR, ...
    matVLST0, matCENTER0, matNPVLST0, matDVE, matADJE, vecROTORBLADES, ...
    valNELE, matROTORHUB, matVEHORIG, vecDVEVEHICLE, vecDVEWING, ...
    matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, vecDVETIP, vecDVELE, ...
    vecDVETE, vecDVEROTORBLADE, vecDVESYM, matROTORAXIS );


matFUSEGEOM = fcnCREATEFUSE(matSECTIONFUSELAGE, vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE);


[ matVEHUVW, matVEHROT, matVEHROTRATE, vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS );
[ matVLST0, matCENTER0, matFUSEGEOM, matROTORHUBGLOB] = fcnROTVEHICLE( matDVE, matVLST0, matCENTER0, valVEHICLES, vecDVEVEHICLE, matVEHORIG, matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS, matROTORHUB, matROTORAXIS, vecROTORVEH);

[ matUINF ] = fcnINITUINF( matCENTER0, matVEHUVW, matVEHROT, vecDVEVEHICLE, ...
    vecDVEROTOR, vecROTORVEH, matVEHORIG, matROTORHUBGLOB, matROTORAXIS, vecROTORRPM );


% update DVE params after vehicle rotation
[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, ~ ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST0);

% [hFig2] = fcnPLOTBODY(0, valNELE, matDVE, matVLST0, matCENTER0,matFUSEGEOM);

%%
% valWSIZE = length(nonzeros(vecDVETE.*(vecDVEWING > 0))); % Amount of wake DVEs shed each timestep
valWSIZE = length(nonzeros(vecDVETE));
%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVESURFACE, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecDVEHVCRD,vecSYM);

%% Alpha Loop

seqALPHA = 0;
seqBETA = 0;
% Preallocating for a turbo-boost in performance
vecCL = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCLF = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCLI = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCDI = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecE = nan(valMAXTIME,valVEHICLES,length(seqALPHA));

for ai = 1:length(seqALPHA)
    
    valALPHA = deg2rad(seqALPHA(ai));
    
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNPVLST = matNPVLST0;
    
    for bi = 1:length(seqBETA)
        
        fprintf('      ANGLE OF ATTACK = %0.3f DEG\n',seqALPHA(ai));
        fprintf('    ANGLE OF SIDESLIP = %0.3f DEG\n',seqBETA(bi));
        fprintf('\n');
        
        valBETA = deg2rad(seqBETA(bi));
        
        % Determining freestream vector
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
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
        vecWDVEWING = [];
        
        % Building wing resultant
        [vecR] = fcnRWING(valNELE, 0, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
            matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            vecWDVETESWP, vecSYM, valWSIZE, flagTRI);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
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
            
            %% Moving the wing
            %             [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE.*(vecDVEWING > 0), matNPVLST);
            %             [matUINF, matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE.*(vecDVEWING > 0), matNPVLST, matFUSEGEOM, vecFUSEVEHICLE);
            [matUINF, matUINFTE, matVEHORIG, ...
                matVLST, matCENTER, ...
                matNEWWAKE, matNPNEWWAKE, ...
                matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, ...
                valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, ...
                vecDVETE, matNPVLST, matFUSEGEOM, vecFUSEVEHICLE, ...
                matVEHROT, vecROTORVEH, matROTORHUBGLOB, ...
                matROTORHUB, matROTORAXIS, vecDVEROTOR, vecROTORRPM );
            %% Generating new wake elements
            %             [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            %                 vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
            %                 = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            %                 vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE.*(vecDVEWING > 0), matWADJE, matNPVLST, vecDVEPANEL, ...
            %                 vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVESURFACE, vecWDVEWING, flagSTEADY, valWSIZE);
            
            [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
                = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, ...
                vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVESURFACE, vecWDVEWING, flagSTEADY, valWSIZE);
            
            if flagPREVIEW ~= 1
                %% Creating and solving WD-Matrix for latest row of wake elements
                % We need to grab from matWADJE only the values we need for this latest row of wake DVEs
                idx = sparse(sum(ismember(matWADJE,[((valWNELE - valWSIZE) + 1):valWNELE]'),2)>0 & (matWADJE(:,2) == 4 | matWADJE(:,2) == 2));
                temp_WADJE = [matWADJE(idx,1) - (valTIMESTEP-1)*valWSIZE matWADJE(idx,2) matWADJE(idx,3) - (valTIMESTEP-1)*valWSIZE];
                
                [matWD, vecWR] = fcnWDWAKE([1:valWSIZE]', temp_WADJE, vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVESYM(end-valWSIZE+1:end), vecWDVETIP(end-valWSIZE+1:end), vecWKGAM(end-valWSIZE+1:end));
                [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, vecWR, valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end-valWSIZE+1:end));
                
                %% Rebuilding and solving wing resultant
                [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
                    matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                    vecWDVETESWP, vecSYM, valWSIZE, flagTRI);
                
                [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
                
                %% Creating and solving WD-Matrix
                [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
                
                %% Relaxing wake
                if valTIMESTEP > 2 && flagRELAX == 1
                    
                    [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
                        vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
                        matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKE(matUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
                        matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVETE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
                        vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
                        vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING);
                    
                    % Creating and solving WD-Matrix
                    [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                    [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
                end
                
                %% Forces
                [vecCL(valTIMESTEP,:,ai), vecCLF(valTIMESTEP,:,ai), vecCLI(valTIMESTEP,:,ai), vecCDI(valTIMESTEP,:,ai), vecE(valTIMESTEP,:,ai), vecDVENFREE, vecDVENIND, ...
                    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, ...
                    vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE, valWNELE, matWDVE, matWVLST, ...
                    matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, ...
                    vecSYM, vecDVETESWP, vecAREA, vecSPAN, [], vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEVEHICLE, valVEHICLES, matVEHROT, flagTRI);
                
            end
            
            %% Post-timestep outputs 
            if flagPRINT == 1 && valTIMESTEP == 1
                
                header1 = ['             '];
                header2 = [' ',sprintf('TIMESTEP'),'    '];
                header3 = ['------------'];
                for j = 1:valVEHICLES
                    header1 = [header1,sprintf('VEHICLE %d',j),'                '];
                    header2 = [header2,[sprintf('CL'),'          ',sprintf('CDI'),'          ']];
                    header3 = [header3,['-------------------------']];
                end
                
                disp(header1);
                disp(header2);
                disp(header3);
            end
            
            txtout = ['  ', sprintf('%4d',valTIMESTEP),'       '];
            for j = 1:valVEHICLES
                txtout = [txtout, sprintf('%0.5f',vecCL(valTIMESTEP,j,ai)), '     ', sprintf('%0.5f',vecCDI(valTIMESTEP,j,ai)), '      '];
            end
            disp(txtout)
            
            
            
            
            if flagGIF == 1
                [hFig3] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM);
                [hFig3] = fcnPLOTWAKE(flagVERBOSE, hFig3, valWNELE, matWDVE, matWVLST, matWCENTER);
                view([33 22])
                
                frame = getframe(hFig3);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                % Write to the GIF File
                
                if valTIMESTEP == 1
                    imwrite(imind,cm,'GIF/output.gif','gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,'GIF/output.gif','gif','WriteMode','append');
                end
            end
            
        end
        
        %% Viscous wrapper
        
        %         [vecCLv(1,ai), vecCD(1,ai), vecPREQ(1,ai), valVINF(1,ai), valLD(1,ai)] = fcnVISCOUS(vecCL(end,ai), vecCDI(end,ai), ...
        %             valWEIGHT, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
        %             vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
        %             matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
        %             valFPWIDTH, valINTERF, vecDVEROLL);
        
    end
end

fprintf('\n');

%% Plotting

if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER, []);
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    
    if flagPLOTWAKEVEL == 1
        try
            quiver3(matWDVEMP(:,1),matWDVEMP(:,2),matWDVEMP(:,3),matWDVEMPIND(:,1),matWDVEMPIND(:,2),matWDVEMPIND(:,3));
        end
    end
    if flagPLOTUINF == 1
        try
            %             quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matVEHUINF(:,1),matVEHUINF(:,2),matVEHUINF(:,3),'g');
            %             quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matROTORUINF(:,1),matROTORUINF(:,2),matROTORUINF(:,3),'c');
            quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matUINF(:,1),matUINF(:,2),matUINF(:,3),'g');
        end
    end
    
    if flagCIRCPLOT == 1
        for i = 1:valNELE
            endpoint_left = (sum(matVLST([matDVE(i,4); matDVE(i,1)],:),1)./2) - matCENTER(i,:);
            endpoint_right = (sum(matVLST([matDVE(i,2); matDVE(i,3)],:),1)./2) - matCENTER(i,:);
            
            tt = fcnGLOBSTAR([endpoint_left; endpoint_right],repmat(vecDVEROLL(i),2,1), repmat(vecDVEPITCH(i),2,1), repmat(vecDVEYAW(i),2,1));
            
            etas = linspace(tt(1,2),tt(2,2))';
            circ = matCOEFF(i,3).*etas.^2 + matCOEFF(i,2).*etas + matCOEFF(i,1);
            
            pt_loc = [linspace(tt(1,1),tt(2,1))' etas circ];
            
            len = size(circ,1);
            circ_glob = fcnSTARGLOB(pt_loc, repmat(vecDVEROLL(i),len,1), repmat(vecDVEPITCH(i),len,1), repmat(vecDVEYAW(i),len,1));
            circ_glob = circ_glob + matCENTER(i,:);
            hold on
            plot3(circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'-m','LineWidth',3)
            hold off
            
        end
    end
    
    
    %     mx = 4;
    %     my = [-8:0.2:8];
    %     mz = [-2:0.2:2];
    %
    %     [X,Y,Z] = meshgrid(mx,my,mz);
    %     fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
    %
    %     q_ind = fcnINDVEL(fpg,valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM,...
    %         valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
    %
    %     % q_ind = s_ind;
    %     hFig20 = figure(20);
    %     clf(20);
    %     % quiver3(fpg(:,1), fpg(:,2), fpg(:,3), w_ind(:,1), w_ind(:,2), w_ind(:,3))
    %     quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
    %     set(gcf,'Renderer','opengl');
    
    
    %     figure(1);
    %     plot(1:valTIMESTEP, eltime)
    %     xlabel('Timestep','FontSize',15)
    %     ylabel('Time per timestep (s)', 'FontSize',15)
    %     box on
    %     grid on
    %     axis tight
    %
    %     figure(3);
    %     plot(1:valTIMESTEP, ttime)
    %     xlabel('Timestep','FontSize',15)
    %     ylabel('Total time (s)', 'FontSize',15)
    %     box on
    %     grid on
    %     axis tight
    
end

% profreport

% hFig23 = figure(23);
% clf(23)
% 
% A = dlmread('size.txt');
% 
% scatter(A(:,1), A(:,2),'kx');
% 
% grid minor
% box on
% axis tight
% xlabel('Number of Input Points','FontSize',15);
% ylabel('Memory Used in VSIND (Gb)','FontSize',15);


%% Viscous wrapper

% whos