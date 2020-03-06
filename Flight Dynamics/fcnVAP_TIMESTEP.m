function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

WAKE.matWAKEGEOM = [];
WAKE.matNPWAKEGEOM = [];
WAKE.vecWDVEHVSPN = [];
WAKE.vecWDVEHVCRD = [];
WAKE.vecWDVEROLL = [];
WAKE.vecWDVEPITCH = [];
WAKE.vecWDVEYAW = [];
WAKE.vecWDVELESWP = [];
WAKE.vecWDVEMCSWP = [];
WAKE.vecWDVETESWP = [];
WAKE.vecWDVEAREA = [];
WAKE.matWDVENORM = [];
WAKE.matWVLST = [];
WAKE.matWDVE = uint32([]);
WAKE.valWNELE = 0;
WAKE.matWCENTER = [];
WAKE.matWCOEFF = [];
WAKE.vecWK = [];
WAKE.matWADJE = uint32([]);
WAKE.vecWDVEPANEL = uint16([]);
WAKE.valLENWADJE = 0;
WAKE.vecWKGAM = [];
WAKE.vecWDVESYM = uint8([]);
WAKE.vecWDVETIP = uint8([]);
WAKE.vecWDVESURFACE = uint8([]);
WAKE.vecWDVETRI = [];
WAKE.vecWPLOTSURF = uint8([]);

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

% Computing structure distributions if data exists
try 
%     [INPU, SURF] = fcnVEHISTRUCT(COND, INPU, SURF, FLAG);
%     [INPU, SURF] = fcnSTRUCTDIST(INPU, SURF, FLAG); 
    FLAG.STRUCTURE = 1; % Create flag if structure data exists
catch
    FLAG.STRUCTURE = 0; 
end

n = 1;
COND.valGUSTTIME = 1;
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
    
    if FLAG.FLIGHTDYN == 1 && valTIMESTEP >= COND.valSTIFFSTEPS
%         [TRIM, VEHI] = fcnSTABDERIV(COND, INPU, SURF, TRIM, VEHI, valTIMESTEP);
        X(valTIMESTEP,1) = -0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*OUTP.vecCD(end);
        Z(valTIMESTEP,1) = -0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*OUTP.vecCL(end);
        M(valTIMESTEP,1) = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC*OUTP.vecVEHCM(end);
        g = 9.81;
        m = COND.vecVEHWEIGHT/g;
        [t,y] = ode45(@(t,y) fcnLONGDYNAMICS(t,y,X(end),Z(end),M(end),m,g,VEHI.vecIYY),[(valTIMESTEP-1)*COND.valDELTIME valTIMESTEP*COND.valDELTIME],[-VEHI.matVEHUVW(1) -VEHI.matVEHUVW(3) TRIM.perturb(valTIMESTEP-1,3) TRIM.perturb(valTIMESTEP-1,4)]);
        TRIM.perturb(valTIMESTEP,:) = y(end,:);
        VEHI.matVEHUVW(1) = -TRIM.perturb(end,1);
        VEHI.matVEHUVW(3) = -TRIM.perturb(end,2);
    else
        TRIM.perturb(valTIMESTEP,1:4) = 0;
    end
    
    % Bend wing if applicable, else move wing normally
    if FLAG.STRUCTURE == 1
        [SURF, INPU, COND, MISC, VISC, OUTP, FLAG, TRIM, VEHI, n] = fcnMOVESTRUCTURE(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, TRIM, valTIMESTEP, n);
        [INPU, SURF, VEHI] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
    else
        [SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);
    end
    
    % Add in gust velocity if applicable
%     if FLAG.GUSTMODE > 0
%         [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old);
%     end
    
    if max(SURF.vecDVEROTOR) > 0 || FLAG.STRUCTURE == 1
        matD = fcnKINCON(matD(1:(size(matD,1)*(2/3)),:), SURF, INPU, FLAG);
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
        tol = 1;
        count = 0;
        
        % Iterating to converge wake and wing coefficients
%         while ~isnan(tol) && tol > 0.005
%             count = 0;
%             delt = 1;
%             while ~isnan(delt) && delt >= 0.005
%                 int_circ1 = fcnINTCIRC2(matPLEX, matCOEFF, matDVEGRID);
%                 vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matCENTER, matKINCON_DVE, matDVECT, vecWDVESYM);
%                 matCOEFF = fcnSOLVED(matD, vecR, valNELE);
%                 [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
%                 matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
% 
%                 % Update wake coefficients
%                 [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE);
%                 matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
%                 int_circ2 = fcnINTCIRC2(matPLEX, matCOEFF, matDVEGRID);
%                 delt = abs((int_circ2 - int_circ1)./int_circ1);
%                 count = count + 1;
%             end
%             vecTSITER(valTIMESTEP) = count;
%             intcirc1 = fcnINTCIRC(SURF);
            [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
            [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
            
%             WAKE.vecWKGAM(end-WAKE.valWSIZE+1:end,1) = [SURF.matCOEFF(SURF.vecDVETE>0,1) + ((WAKE.vecWDVEHVSPN(end-WAKE.valWSIZE+1:end).^2)./3).*SURF.matCOEFF(SURF.vecDVETE>0,3)];

            %% Creating and solving WD-Matrix
            [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
            [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
            
%             intcirc2 = fcnINTCIRC(SURF);
%             
%             tol = abs((intcirc2-intcirc1)./intcirc1);
%             count = count + 1;
%         end
%         OUTP.vecTSITER(valTIMESTEP,1) = count;
        
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
%         fcnPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.vecCT, MISC.vecROTORJ, VEHI.vecROTORVEH, 1)
        fcnFLIGHTDYNPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.matDEFGLOB, INPU.vecVEHCG, TRIM.perturb, 1)
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

end
