function [FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT_FLEX(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP)

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
    if FLAG.OPT == 1
        [INPU, SURF, VEHI] = fcnVEHISTRUCT_OPT(COND, INPU, SURF, FLAG, VEHI);
        FLAG.STRUCTURE = 1; % Create flag if structure data exists
    elseif FLAG.OPT == 2
        [INPU, SURF, VEHI] = fcnVEHISTRUCT_PARAM(COND, INPU, SURF, FLAG, VEHI);
        FLAG.STRUCTURE = 1; % Create flag if structure data exists
    else
        [INPU, SURF, VEHI] = fcnVEHISTRUCT(COND, INPU, SURF, FLAG, VEHI);
        FLAG.STRUCTURE = 1; % Create flag if structure data exists
    end
catch
    FLAG.STRUCTURE = 0; 
end

SURF.valTBOOM = SURF.matVLST(SURF.matDVE(SURF.idxTAIL(1),1),1) - SURF.matVLST(SURF.matDVE(SURF.idxFLEX(1),1),1);
[INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle

n = 1;
COND.valGUSTTIME = 1;
SURF.gust_vel_old = zeros(SURF.valNELE,1);

end