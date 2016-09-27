clc
clear

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

strFILE = 'VAP input.txt';

[flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF] = fcnVAPREAD(strFILE);

% strFILE = 'input.txt';
%
% [flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF] = fcnFWREAD(strFILE);

flagPLOT = 1;

%% Discretize geometry into DVEs

[matCENTER, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, vecDVENORM, ...
    matVLST, matDVE, valNELE, matADJE, vecDVESYM, vecDVETIP] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);


%% Plotting Wing

if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(1, valNELE, matDVE, matVLST, matCENTER, vecDVENORM);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
end

%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

%% Create D-Resultant, solve D-Matrix

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

%% Viscous wrapper

toc

% whos