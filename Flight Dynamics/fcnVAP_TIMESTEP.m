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
timestep_folder = 'timestep_data\';

%% Add boundary conditions to D-Matrix
[matD] = fcnDWING(SURF, INPU);

%% Add kinematic conditions to D-Matrix
[SURF.vecK] = fcnSINGFCT(SURF.valNELE, SURF.vecDVESURFACE, SURF.vecDVETIP, SURF.vecDVEHVSPN);
[matD, SURF.matCOLLPTS] = fcnKINCON(matD, SURF, INPU, FLAG);

%% Preparing to timestep
% Building wing resultant
[vecR,w_wake] = fcnRWING(0, SURF, WAKE, FLAG);

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

VEHI.matVEHUVW = [COND.vecVEHVINF*cosd(COND.vecVEHALPHA),0,COND.vecVEHVINF*sind(COND.vecVEHALPHA)];
SURF.GammaInt = [];

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
    SURF.matCENTER_t(:,:,valTIMESTEP) = SURF.matCENTER;
    SURF.matUINFdir = -SURF.matUINF./sqrt(sum(SURF.matUINF.^2,2));
    UINFdir = repmat(VEHI.matGLOBUVW./vecnorm(VEHI.matGLOBUVW,2,2),size(SURF.matUINFdir,1),1);
    for i = 1:size(SURF.matUINFdir,1)
        VEHI.ddir(i,:) = ([-1 0 0; 0 1 0; 0 0 -1]*SURF.matUINFdir(i,:)')'; % Drag direction is parallel to Uinf
        VEHI.ldir(i,:) = ([0 0 -1; 0 1 0; 1 0 0]*VEHI.ddir(i,:)')'; % Lift direction is perpendicular to Uinf
        ddir(i,:) = ([-1 0 0; 0 1 0; 0 0 -1]*UINFdir(i,:)')'; % Drag direction is parallel to Uinf
        ldir(i,:) = ([0 0 -1; 0 1 0; 1 0 0]*VEHI.ddir(i,:)')'; % Lift direction is perpendicular to Uinf
    end
    
    OUTP.matGLOBUVW(valTIMESTEP,:) = VEHI.matGLOBUVW;
        
    if FLAG.FLIGHTDYN == 1 && valTIMESTEP > 2
        
        iter = 0;
        g = 9.81;
        
        % Determine aerodynamic force directions and vehicle flight path
        % AoA and flight path angle
        COND.vecVEHALPHA = atand(VEHI.matVEHUVW(3)/VEHI.matVEHUVW(1));
        OUTP.vecVEHALPHA(valTIMESTEP,1) = COND.vecVEHALPHA;
        COND.vecVEHFPA = -atand(VEHI.matGLOBUVW(3)/VEHI.matGLOBUVW(1));
        OUTP.vecVEHFPA(valTIMESTEP,1) = COND.vecVEHFPA;
%         COND.valDELTIME = INPU.vecCMAC/INPU.vecM(1)/COND.vecVEHVINF;
        if valTIMESTEP > 4
%             OUTP.GlobForce(valTIMESTEP-1,:) = OUTP.GlobForce(valTIMESTEP-1,:) - m*SURF.matBEAMACC(valTIMESTEP-1,:);
        end
        
        if FLAG.STIFFWING == 0
            
            [OUTP.BForceFuse(valTIMESTEP,:)] = fcnGLOBSTAR(OUTP.GlobForceFuse(valTIMESTEP-1,:), 0, pi+VEHI.vecVEHDYN(valTIMESTEP-1,4), 0); % Rotation of force in earth frame to body frame
            OUTP.BForceFuse(valTIMESTEP,3) = OUTP.BForceFuse(valTIMESTEP,3) - 2*OUTP.vecFUSEV(valTIMESTEP-1,1);
            
            [OUTP.BForce(valTIMESTEP,:)] = fcnGLOBSTAR(OUTP.GlobForce(valTIMESTEP-1,:), 0, pi+VEHI.vecVEHDYN(valTIMESTEP-1,4), 0); % Rotation of force in earth frame to body frame
            
            Mstruct = 2*OUTP.vecFUSEM(valTIMESTEP-1,1) + 2*OUTP.vecFUSEV(valTIMESTEP-1,1)*(INPU.vecVEHCG(1) - SURF.matBEAMLOC(1,1,valTIMESTEP-1));
            [OUTP] = fcnTAILMOM(SURF, VEHI, OUTP, INPU, COND, FLAG);
            M = Mstruct*COND.valDENSITY + OUTP.vecVEHPITCHMOM;
%             M = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC*OUTP.vecVEHCM;
            
            m = COND.vecVEHWEIGHT/g;
            mf = VEHI.vecWINGMASS(2) + sum(VEHI.vecFUSEMASS,1) + VEHI.vecPAYLMASS;
            
            [t,y] = ode45(@(t,y) fcnFUSEDYNAMICS(t,y,OUTP.BForce(valTIMESTEP,1),OUTP.BForceFuse(valTIMESTEP,3),M,m,mf,g,VEHI.vecIYY),...
                [0 COND.valDELTIME],[VEHI.matVEHUVW(1) VEHI.matVEHUVW(3) VEHI.vecVEHDYN(valTIMESTEP-1,3) VEHI.vecVEHDYN(valTIMESTEP-1,4)]);
            VEHI.vecVEHDYN(valTIMESTEP,:) = y(end,:);
        end
                         
%         [OUTP.BForce(valTIMESTEP,:)] = fcnGLOBSTAR(OUTP.GlobForce(valTIMESTEP-1,:), 0, pi+VEHI.vecVEHDYN(valTIMESTEP-1,4), 0); % Rotation of force in earth frame to body frame
%         
%         M(valTIMESTEP,1) = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA*INPU.vecCMAC*OUTP.vecVEHCM;
%         
%         m = COND.vecVEHWEIGHT/g;
%         mf = VEHI.vecWINGMASS(2) + VEHI.vecFUSEMASS + VEHI.vecPAYLMASS;
%         [t,y] = ode23t(@(t,y) fcnLONGDYNAMICS(t,y,OUTP.BForce(valTIMESTEP,1),OUTP.BForce(valTIMESTEP,3),M(end),m,mf,g,VEHI.vecIYY),[0 COND.valDELTIME],[VEHI.matVEHUVW(1) VEHI.matVEHUVW(3) VEHI.vecVEHDYN(valTIMESTEP-1,3) VEHI.vecVEHDYN(valTIMESTEP-1,4)]);
%         VEHI.vecVEHDYN(valTIMESTEP,:) = y(end,:);
              
        % Vehicle velocity in the body-fixed frame
        VEHI.matVEHUVW(1) = VEHI.vecVEHDYN(valTIMESTEP,1);
        VEHI.matVEHUVW(3) = VEHI.vecVEHDYN(valTIMESTEP,2);
        
        [VEHI.matGLOBUVW] = fcnSTARGLOB(VEHI.matVEHUVW, 0, pi+VEHI.vecVEHDYN(valTIMESTEP,4), 0); % Rotation of velocity in body frame to earth frame
        
        OUTP.matGLOBUVW(valTIMESTEP,:) = VEHI.matGLOBUVW;
        
        OUTP.vecVEHORIG(valTIMESTEP,1) = INPU.matVEHORIG(3);
        
        if valTIMESTEP > 3
            SURF.matBEAMACC(valTIMESTEP,:) = (OUTP.matGLOBUVW(valTIMESTEP,:) - OUTP.matGLOBUVW(valTIMESTEP-1,:))./COND.valDELTIME;
        end
        
        COND.vecVEHVINF = sqrt(sum(VEHI.matVEHUVW.^2)); % Vehicle Uinf magnitude
        OUTP.vecVEHVINF(valTIMESTEP,1) = COND.vecVEHVINF;
    else
        VEHI.vecVEHDYN(valTIMESTEP,1) = VEHI.matVEHUVW(1);
        VEHI.vecVEHDYN(valTIMESTEP,2) = VEHI.matVEHUVW(3);
        VEHI.vecVEHDYN(valTIMESTEP,3) = 0;
        VEHI.vecVEHDYN(valTIMESTEP,4) = deg2rad(COND.vecVEHPITCH);
    end
    
    OUTP.dt(valTIMESTEP,1) = COND.valDELTIME;
    OUTP.sim_time = cumsum(OUTP.dt,1); % Store updated simulation time based on new dt
    
    % Bend wing if applicable, else move wing normally
    if FLAG.STRUCTURE == 1
        [INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
        [SURF, INPU, COND, MISC, VISC, OUTP, FLAG, TRIM, VEHI, n] = fcnMOVESTRUCTURE(INPU, VEHI, MISC, COND, SURF, VISC, FLAG, OUTP, TRIM, valTIMESTEP, n);
        [INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND); % Recompute mass properties of vehicle
        OUTP.vecCGLOC(valTIMESTEP,:) = INPU.vecVEHCG;
        
        tempCG = fcnGLOBSTAR(OUTP.vecCGLOC(valTIMESTEP,:),0,VEHI.vecVEHDYN(valTIMESTEP,4),0);
        tempORIG = fcnGLOBSTAR(INPU.matVEHORIG,0,VEHI.vecVEHDYN(valTIMESTEP,4),0);
        OUTP.vecCGLOC_pC(valTIMESTEP,:) = tempCG(1) - tempORIG(1);
        VEHI.vecFUSECG_t(valTIMESTEP,:) = VEHI.vecFUSECG;
        OUTP.vecWINGCG(valTIMESTEP,:) = VEHI.vecWINGCG(1,:);
        
        % Specific kinetic, potential and total vehicle energy
        OUTP.vecVEHENERGY(valTIMESTEP,1) = 0.5*COND.vecVEHVINF*COND.vecVEHVINF;
        OUTP.vecVEHENERGY(valTIMESTEP,2) = 9.81*(OUTP.vecCGLOC(valTIMESTEP,3)-OUTP.vecCGLOC(1,3));
        OUTP.vecVEHENERGY(valTIMESTEP,3) = OUTP.vecVEHENERGY(valTIMESTEP,1) + OUTP.vecVEHENERGY(valTIMESTEP,2);
        
        if (FLAG.STIFFWING == 0 && FLAG.FLIGHTDYN == 1) || FLAG.FLIGHTDYN == 1
            % Calculate vehicle energy state
            [OUTP] = fcnVEHENERGY(INPU, COND, SURF, OUTP, VEHI, FLAG, valTIMESTEP);
            if valTIMESTEP >= COND.valGUSTSTART
                OUTP.vecZE_old(valTIMESTEP,1) = (OUTP.vecVEHVINF(end)^2 - OUTP.vecVEHVINF(COND.valGUSTSTART)^2)/(2*9.81) + OUTP.vecCGLOC(end,3) - OUTP.vecCGLOC(COND.valGUSTSTART,3);
            end
        end
        
        OUTP.vecTIPPITCH(valTIMESTEP,1) = SURF.vecDVEPITCH(INPU.vecN(1));
    else
        [SURF, INPU, MISC, VISC] = fcnMOVESURFACE(INPU, VEHI, MISC, COND, SURF, VISC);
    end
    
    % Add in gust velocity if applicable
%     if FLAG.GUSTMODE > 0
%         [SURF.matUINF, SURF.gust_vel_old] = fcnGUSTWING(SURF.matUINF,COND.valGUSTAMP,COND.valGUSTL,FLAG.GUSTMODE,COND.valDELTIME,COND.vecVEHVINF,COND.valGUSTSTART,SURF.matCENTER,SURF.gust_vel_old,COND.start_loc);
%     end   
    
    if max(SURF.vecDVEROTOR) > 0 || FLAG.STRUCTURE == 1
        [matD, SURF.matCOLLPTS] = fcnKINCON(matD(1:(size(matD,1)*(2/3)),:), SURF, INPU, FLAG);
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
        
                        
        if valTIMESTEP == 1
            OUTP.DEBUG.vecWKGAM = [];
            OUTP.DEBUG.w_wake = [];
        end
        
        if valTIMESTEP > 1
            OUTP.DEBUG.vecWKGAM(size(OUTP.DEBUG.vecWKGAM,1)+1:size(WAKE.vecWKGAM,1),valTIMESTEP-1) = 0;
        end
        OUTP.DEBUG.vecWKGAM(:,valTIMESTEP) = WAKE.vecWKGAM;
        OUTP.DEBUG.vecR(:,valTIMESTEP) = vecR;
        OUTP.DEBUG.matDVENORM(:,:,valTIMESTEP) = SURF.matDVENORM;
        OUTP.DEBUG.w_wake(:,:,valTIMESTEP) = w_wake;
        OUTP.DEBUG.matUINF(:,:,valTIMESTEP) = SURF.matUINF;
        
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
                
                [vecR,w_wake] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
                [SURF.matCOEFF] = fcnSOLVED(matD, vecR, SURF.valNELE);
                
                WAKE.vecWKGAM = fcnWKGAM(FLAG.STEADY, WAKE.vecWKGAM, SURF.matCOEFF, SURF.vecDVETE, SURF.vecDVEHVSPN(SURF.vecDVETE > 0), WAKE.valWNELE, WAKE.valWSIZE);
                
                [matWD, WAKE.vecWR] = fcnWDWAKE([1:WAKE.valWNELE]', WAKE.matWADJE, WAKE.vecWDVEHVSPN, WAKE.vecWDVESYM, WAKE.vecWDVETIP, WAKE.vecWKGAM, INPU.vecN);
                [WAKE.matWCOEFF] = fcnSOLVEWD(matWD, WAKE.vecWR, WAKE.valWNELE, WAKE.vecWKGAM, WAKE.vecWDVEHVSPN);
                
                delt = max(max(abs(SURF.matCOEFF - matCOEFF1)));
                
                if valTIMESTEP > 1
                    OUTP.DEBUG.vecWKGAM(size(OUTP.DEBUG.vecWKGAM,1)+1:size(WAKE.vecWKGAM,1),valTIMESTEP-1) = 0;
                end
                OUTP.DEBUG.vecWKGAM(:,valTIMESTEP) = WAKE.vecWKGAM;
                OUTP.DEBUG.vecR(:,valTIMESTEP) = vecR;
                OUTP.DEBUG.matDVENORM(:,:,valTIMESTEP) = SURF.matDVENORM;
                OUTP.DEBUG.w_wake(:,:,valTIMESTEP) = w_wake;
                OUTP.DEBUG.matUINF(:,:,valTIMESTEP) = SURF.matUINF;
%                 disp(delt)
            end
        else
            [vecR,w_wake] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
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
                    
                    [vecR,~] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG);
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
            OUTP.GlobForce(valTIMESTEP,:) = 2*COND.valDENSITY*sum(dot(SURF.matDVEIFORCE,VEHI.ldir,2).*VEHI.ldir + SURF.matDVEINDDRAG.*VEHI.ddir,1);
            OUTP.GlobForceFuse(valTIMESTEP,:) = 2*COND.valDENSITY*(sum(dot(SURF.matDVEIFORCE(SURF.vecDVEWING == 2,:),VEHI.ldir(SURF.vecDVEWING == 2,:),2).*VEHI.ldir(SURF.vecDVEWING == 2,:),1)...
                + sum(SURF.matDVEINDDRAG.*VEHI.ddir,1));
            OUTP.matUINF_t(:,:,valTIMESTEP) = SURF.matUINF;
%             SURF.matBEAMLOC(:,:,valTIMESTEP) = SURF.matEALST(1:INPU.valNSELE,:);
        end
        
        if FLAG.SAVETIMESTEP == 1
            save([timestep_folder, 'timestep_', num2str(valTIMESTEP), '.mat'],'valTIMESTEP','INPU','COND','MISC','WAKE','VEHI','SURF','OUTP','VISC','TRIM','FLAG','iter');
        end
    end
    
    %% Post-timestep outputs
    if FLAG.PRINT == 1
%         fcnPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.vecCT, MISC.vecROTORJ, VEHI.vecROTORVEH, 1)
        fcnFLIGHTDYNPRINTOUT(FLAG.PRINT, valTIMESTEP, INPU.valVEHICLES, OUTP.vecCL, OUTP.vecCDI, OUTP.matDEFGLOB, INPU.vecVEHCG, VEHI.vecVEHDYN, 1)
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
