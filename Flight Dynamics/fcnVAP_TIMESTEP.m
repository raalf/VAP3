function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC, iter)

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
    
    if FLAG.FLIGHTDYN == 1 && valTIMESTEP > COND.valSTIFFSTEPS
        
        SURF.matUINFdir = SURF.matUINF./sqrt(sum(abs(SURF.matUINF.^2),2));
%         SURF.matUINFdir = VEHI.matVEHUVW./sqrt(sum(abs(VEHI.matVEHUVW.^2),2));
        for i = 1:size(SURF.matUINFdir,1)
            ddir(i,:) = ([1 0 0; 0 1 0; 0 0 1]*SURF.matUINFdir(i,:)')';
            ldir(i,:) = ([0 0 -1; 0 1 0; 1 0 0]*ddir(i,:)')';
        end
        lift = COND.valDENSITY*(SURF.vecDVELFREE+SURF.vecDVELIND).*ldir;
        drag = COND.valDENSITY*SURF.vecDVEDIND.*ddir;
        GlobForce = 2*sum(lift + drag,1) - 0*0.5*COND.valDENSITY*10*10*INPU.vecAREA*TRIM.valCDI.*ddir(1,:);

        COND.vecVEHVINF = sqrt(sum(VEHI.matVEHUVW.^2));
        COND.valDELTIME = 0.5/COND.vecVEHVINF;
        OUTP.dt(valTIMESTEP,1) = COND.valDELTIME;
        OUTP.sim_time = cumsum(OUTP.dt,1); % Store updated simulation time based on new dt
%         vdir = VEHI.matGLOBUVW./sqrt(sum(abs(VEHI.matGLOBUVW.^2))); % Global velocity direction vector
%         ddir = ([-1 0 0; 0 1 0; 0 0 -1]*vdir')'; % Drag force direction vector
%         ldir = ([0 0 -1; 0 1 0; 1 0 0]*ddir')'; % Lift force direction vector
%         drag = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*(OUTP.vecCDI(valTIMESTEP-1)).*ddir - 0.5*COND.valDENSITY*10*10*INPU.vecAREA*TRIM.valCDI.*ddir; % Drag force vector in global frame
%         lift = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*OUTP.vecCL(valTIMESTEP-1)*ldir; % Lift force vector in global frame
        
        VEHI.matROTMAT = [cos(pi-TRIM.perturb(valTIMESTEP-1,4)) 0 -sin(pi-TRIM.perturb(valTIMESTEP-1,4));...
                  0   1   0;...
                  sin(pi-TRIM.perturb(valTIMESTEP-1,4)) 0 cos(pi-TRIM.perturb(valTIMESTEP-1,4))]; % Rotation matrix from body frame to earth frame
        
%         BForce(valTIMESTEP,:) = (VEHI.matROTMAT'*drag' + VEHI.matROTMAT'*lift')'; % Aerodynamic force vector in the vehicle body frame
        OUTP.BForce(valTIMESTEP,:) = (VEHI.matROTMAT'*GlobForce')'; % Aerodynamic force vector in the vehicle body frame
%         OUTP.BForce(valTIMESTEP,:) = (VEHI.matROTMAT*GlobForce')'; % Aerodynamic force vector in the vehicle body frame
%         OUTP.BForce(valTIMESTEP,:) = GlobForce; % Aerodynamic force vector in the vehicle body frame
                              
        M(valTIMESTEP,1) = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC*OUTP.vecVEHCM;
        g = 9.81;
        m = COND.vecVEHWEIGHT/g;
        [t,y] = ode45(@(t,y) fcnLONGDYNAMICS(t,y,OUTP.BForce(valTIMESTEP,1),OUTP.BForce(valTIMESTEP,3),M(end),m,g,VEHI.vecIYY),[0 COND.valDELTIME],[VEHI.matVEHUVW(1) VEHI.matVEHUVW(3) TRIM.perturb(valTIMESTEP-1,3) TRIM.perturb(valTIMESTEP-1,4)]);
        TRIM.perturb(valTIMESTEP,:) = y(end,:);
        
        % Vehicle velocity in the body-fixed frame
        VEHI.matVEHUVW(1) = TRIM.perturb(valTIMESTEP,1);
        VEHI.matVEHUVW(3) = TRIM.perturb(valTIMESTEP,2);
        
        VEHI.matROTMAT = [cos(pi-TRIM.perturb(valTIMESTEP,4)) 0 -sin(pi-TRIM.perturb(valTIMESTEP,4));...
                          0   1   0;...
                          sin(pi-TRIM.perturb(valTIMESTEP,4)) 0 -cos(pi-TRIM.perturb(valTIMESTEP,4))]; % Rotation matrix from body frame to earth frame
        VEHI.matGLOBUVW = (VEHI.matROTMAT*VEHI.matVEHUVW')'; % Vehicle velocity in the earth frame
        OUTP.matGLOBUVW(valTIMESTEP,:) = VEHI.matGLOBUVW;
    else
        TRIM.perturb(valTIMESTEP,1:4) = 0;
    end
    
    % Bend wing if applicable, else move wing normally
    if FLAG.STRUCTURE == 1
        [INPU, SURF, VEHI] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
        [SURF, INPU, COND, MISC, VISC, OUTP, FLAG, TRIM, VEHI, n] = fcnMOVESTRUCTURE(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, TRIM, valTIMESTEP, n);
        [INPU, SURF, VEHI] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
        OUTP.vecCGLOC(valTIMESTEP,:) = INPU.vecVEHCG;
        
        % Specific kinetic, potential and total vehicle energy
        OUTP.vecVEHENERGY(valTIMESTEP,1) = 0.5*COND.vecVEHVINF*COND.vecVEHVINF;
        OUTP.vecVEHENERGY(valTIMESTEP,2) = 9.81*(OUTP.vecCGLOC(valTIMESTEP,3)-OUTP.vecCGLOC(1,3));
        OUTP.vecVEHENERGY(valTIMESTEP,3) = OUTP.vecVEHENERGY(valTIMESTEP,1) + OUTP.vecVEHENERGY(valTIMESTEP,2);
    else
        [SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);
    end
    
    % Add in gust velocity if applicable
    if FLAG.GUSTMODE > 0
        [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old);
    end   
    
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
        
%         %% Rebuilding and solving wing resultant        
%         [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
%         [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
% 
%         %% Creating and solving WD-Matrix
%         [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
%         [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
            
        %% Iteration Block #1
        if iter == true
            delt = 1;
            while delt >= 0.000001
                matCOEFF1 = SURF.matCOEFF;
                
                [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
                [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
                
                WAKE.vecWKGAM = fcnWKGAM(FLAG.STEADY, WAKE.vecWKGAM, SURF.matCOEFF, SURF.vecDVETE, SURF.vecDVEHVSPN(SURF.vecDVETE > 0), WAKE.valWNELE, WAKE.valWSIZE);
                
                [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
                [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
                
                delt = max(max(abs(SURF.matCOEFF - matCOEFF1)));
%                 disp(delt)
            end
        else
            [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
            [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
            
            [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
            [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
        end
        
        %% Relaxing wake
        if valTIMESTEP > 2 && FLAG.RELAX == 1
            old_span = WAKE.vecWDVEHVSPN;
            WAKE = fcnRELAXWAKE(valTIMESTEP, SURF, WAKE, COND, FLAG, INPU);
            WAKE.matWCOEFF(:,2:3) = WAKE.matWCOEFF(:,2:3).*[old_span./WAKE.vecWDVEHVSPN (old_span./WAKE.vecWDVEHVSPN).^2];
            
            %% Iteration Block #2
            if iter == true
                delt = 1;
                while delt >= 0.000001
                    matCOEFF1 = SURF.matCOEFF;
                    
                    [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
                    [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
                    
                    WAKE.vecWKGAM = fcnWKGAM(FLAG.STEADY, WAKE.vecWKGAM, SURF.matCOEFF, SURF.vecDVETE, SURF.vecDVEHVSPN(SURF.vecDVETE > 0), WAKE.valWNELE, WAKE.valWSIZE);
                    
                    [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
                    [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
                    
                    delt = max(max(abs(SURF.matCOEFF - matCOEFF1)));
%                     disp(delt)
                end
            end
            
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
