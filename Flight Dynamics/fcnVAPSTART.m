function [FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP] = fcnVAPSTART(filename,VAP_IN)
% Function that performs all the initialization tasks such as reading in
% data and generating DVEs

warning off
if nargin == 0
    VAP_MAIN;
    return
end

%% Reading in geometry
[FLAG, COND, VISC, INPU, VEHI, SURF] = fcnXMLREAD(filename, VAP_IN);

FLAG.GPU = 0;
FLAG.PRINT = 1;
FLAG.PLOT = 1;
FLAG.VISCOUS = 0;
FLAG.CIRCPLOT = 0;
FLAG.GIF = 0;
FLAG.PREVIEW = 0;
FLAG.PLOTWAKEVEL = 0;
FLAG.PLOTUINF = 0;
FLAG.VERBOSE = 0;
FLAG.SAVETIMESTEP = 0;
FLAG.TRIMMED = 0;
FLAG.STATICAERO = 0;
FLAG.NACELLE = 0;

% Initializing parameters to null/zero/nan
[WAKE, OUTP, INPU, SURF] = fcnINITIALIZE(COND, INPU, SURF);

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
    if exist('timesteps/') ~= 7; mkdir(timesteps); end
    timestep_folder = ['timesteps/',regexprep(filename,{'inputs/', '.vap'}, ''), '_(', datestr(now, 'dd_mm_yyyy HH_MM_SS_FFF'),')/'];
    mkdir(timestep_folder);
end

% Check if the files required by the viscous calculations exist
[FLAG] = fcnVISCOUSFILECHECK(FLAG, VISC);

% %% Discretizing geometry into DVEs
% % Adding collective pitch to the propeller/rotor
% if ~isempty(COND.vecCOLLECTIVE)
%     INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) = INPU.matGEOM(:,5,INPU.vecPANELROTOR > 0) + repmat(reshape(COND.vecCOLLECTIVE(INPU.vecPANELROTOR(INPU.vecPANELROTOR > 0), 1),1,1,[]),2,1,1);
% end
% [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, FLAG, OUTP, SURF);
% % fcnPLOTPKG([], FLAG, SURF, VISC, WAKE, COND, INPU)
% %% Advance Ratio
% MISC.vecROTORJ = [];
% for jj = 1:length(COND.vecROTORRPM)
%     MISC.vecROTORJ(jj) = (COND.vecVEHVINF(VEHI.vecROTORVEH(jj))*60)./(abs(COND.vecROTORRPM(jj)).*INPU.vecROTDIAM(jj));
% end
% 
% %% Add boundary conditions to D-Matrix
% [matD] = fcnDWING(SURF, INPU);
% 
% %% Add kinematic conditions to D-Matrix
% [SURF.vecK] = fcnSINGFCT(SURF.valNELE, SURF.vecDVESURFACE, SURF.vecDVETIP, SURF.vecDVEHVSPN);
% [matD] = fcnKINCON(matD, SURF, INPU, FLAG);
% 
% %% Preparing to timestep
% % Building wing resultant
% [vecR] = fcnRWING(0, SURF, WAKE, FLAG);
% 
% % Solving for wing coefficients
% [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
% 
% SURF.matNPDVE = SURF.matDVE;
% % Computing structure distributions if data exists
% try 
%     [INPU, SURF] = fcnVEHISTRUCT(INPU, SURF, FLAG);
% %     [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF, FLAG); 
%     FLAG.STRUCTURE = 1; % Create flag if structure data exists
% catch
%     FLAG.STRUCTURE = 0; 
% end
% 
% n = 1;
% COND.valGUSTTIME = 1;
% SURF.gust_vel_old = zeros(SURF.valNELE,1);


end