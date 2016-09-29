clc
clear
clf
tic

disp('===========================================================================');
disp('VAP (Based on FreeWake 2015)');
disp('Running Version 2016.09  .                             .');
disp('Includes stall model    //                             \\');
disp('No trim solution       //                               \\');
disp('                      //                                 \\');
disp('                     //                _._                \\');
disp('                  .---.              .//|\\.              .---.');
disp('         ________/ .-. \_________..-~ _.-._ ~-..________ / .-. \_________');
disp('                 \ ~-~ /   /H-     `-=.___.=-''     -H\   \ ~-~ /');
disp('                   ~~~    / H          [H]          H \    ~~~');
disp('                         / _H_         _H_         _H_ \');
disp('                           UUU         UUU         UUU');
disp('===========================================================================');
disp(' ');

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Reading in geometry

% strFILE = 'VAP christmas.txt';
% % strFILE = 'VAP input.txt';
% [flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF] = fcnVAPREAD(strFILE);

strFILE = 'input.txt';

[flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF] = fcnFWREAD(strFILE);

flagPLOT = 1;

%% Discretize geometry into DVEs

[matCENTER, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVETE] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);


%% Plotting Wing

if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(1, valNELE, matDVE, matVLST, matCENTER, matDVENORM);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
end

%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER, matVLST, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM);


%% Alpha Loop
for ai = 1:length(seqALPHA)
    valALPHA = deg2rad(seqALPHA(ai));
    for bi = 1:length(seqBETA)
        valBETA = deg2rad(seqBETA(bi));
        
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
        % Building wing resultant
        [vecR] = fcnRWING(valNELE, 0, matCENTER, matDVENORM, vecUINF);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
        matWAKEGEOM = [];
        for valTIMESTEP = 1:1%valMAXTIME
            %% Timestep to solution
            %   Move wing
            %   Generate new wake elements
            %   Create W-Matrix and W-Resultant
            %   Solve W-Matrix
            %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
            %   Calculate surface normal forces
            %   Calculate DVE normal forces
            %   Calculate induced drag
            %   Calculate cn, cl, cy, cdi
            %   Calculate viscous effects
            
            % Moving the wing
            [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE);
            
            
            
            matWCENTER = mean(matNEWWAKE,3);
            [ vecWDVEHVSPN, vecWDVEHVCRD, ...
                vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
                vecWDVELESWP, vecWDVEMCSWP, vecWDVETESWP, ...
                vecWDVEAREA, matWDVENORM, ...
                matWVLST, matWDVE, valWNELE ] = fcnDVECORNER2PARAM(matWCENTER, matNEWWAKE(:,:,1), matNEWWAKE(:,:,2), matNEWWAKE(:,:,3), matNEWWAKE(:,:,4));
            clf(2)
            [hFig2] = fcnPLOTBODY(1, valNELE, matDVE, matVLST, matCENTER, matDVENORM);
            [hFig2] = fcnPLOTBODY(0, 0, matWDVE, matWVLST, matWCENTER, matWDVENORM)
            % Generating new wake elements
            %             [matWAKEGEOM, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM);
            
        end
    end
end

%% Viscous wrapper

toc

% whos