clc
clear

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

strFILE = 'inputs/VAP christmas.txt';
% strFILE = 'inputs/VAP input.txt';

[flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF] = fcnVAPREAD(strFILE);

vecM = [3 3 3 3 3 3 2]';
seqALPHA = [5 7];

% strFILE = 'inputs/input.txt';
% strFILE = 'inputs/Config 1.txt';
% strFILE = 'inputs/Config 2.txt';

% [flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF] = fcnFWREAD(strFILE);


valMAXTIME = 10;
flagRELAX = 1;

flagPLOT = 1;
flagVERBOSE = 1;

%% Discretize geometry into DVEs

[matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, ...
    vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

valWSIZE = length(nonzeros(vecDVETE)); % Amount of wake DVEs shed each timestep

%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecSYM);

%% Alpha Loop

% Preallocating for a turbo-boost in performance
vecCL = zeros(valMAXTIME, length(seqALPHA));
vecCDI = zeros(valMAXTIME, length(seqALPHA));
vecE = zeros(valMAXTIME, length(seqALPHA));

for ai = 1:length(seqALPHA);
    fprintf('\nAlpha = %0.3f', seqALPHA(ai));
    valALPHA = deg2rad(seqALPHA(ai));
    
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNPVLST = matNPVLST0;
    
    for bi = 1:length(seqBETA)
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
        [vecR] = fcnRWING(valNELE, 0, matCENTER, matDVENORM, vecUINF, valWNELE, matWDVE, ...
            matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            vecWDVETESWP, vecSYM, valWSIZE);
        
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
            [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNPVLST);
            
            %% Generating new wake elements
            [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
                = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, ...
                vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVEWING, vecWDVEWING, flagSTEADY, valWSIZE);
            
            %% Creating and solving WD-Matrix for latest row of wake elements
            % We need to grab from matWADJE only the values we need for this latest row of wake DVEs
            idx = sparse(sum(ismember(matWADJE,[((valWNELE - valWSIZE) + 1):valWNELE]'),2)>0 & (matWADJE(:,2) == 4 | matWADJE(:,2) == 2));
            temp_WADJE = [matWADJE(idx,1) - (valTIMESTEP-1)*valWSIZE matWADJE(idx,2) matWADJE(idx,3) - (valTIMESTEP-1)*valWSIZE];
            
            [matWD, vecWR] = fcnWDWAKE([1:valWSIZE]', temp_WADJE, vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVESYM(end-valWSIZE+1:end), vecWDVETIP(end-valWSIZE+1:end), vecWKGAM(end-valWSIZE+1:end));
            [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, vecWR, valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end-valWSIZE+1:end));
            
            %% Rebuilding and solving wing resultant
            [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, vecUINF, valWNELE, matWDVE, ...
                matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVETESWP, vecSYM, valWSIZE);
            
            [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
            
            %% Creating and solving WD-Matrix
            [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
            [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
            
            %% Relaxing wake
            if valTIMESTEP > 2 && flagRELAX == 1
                
                [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
                    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWDVENORM, ...
                    matWVLST, matWDVE, idxWVLST, vecWK] = fcnRELAXWAKE(matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
                    matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVEHVSPN, vecDVELESWP, ...
                    vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVELESWP, vecWDVEPITCH, ...
                    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING);
                
                % Creating and solving WD-Matrix
                [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
            end
            
            %% Timing
            %             eltime(valTIMESTEP) = toc;
            %             ttime(valTIMESTEP) = sum(eltime);
            
            %% Forces
            
            [vecCL(valTIMESTEP,ai), vecCDI(valTIMESTEP,ai), vecE(valTIMESTEP,ai), vecDVENFREE, vecDVENIND, ...
                vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = ...
                fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, vecUINF, vecDVELESWP,...
                vecDVEMCSWP, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE,...
                valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, ...
                vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP, valAREA, valSPAN, valBETA, ...
                vecDVEWING, vecN, vecM, vecDVEPANEL);
            
%             fprintf('\n\tTimestep = %0.0f', valTIMESTEP);
%             fprintf('\tCL = %0.5f',vecCL(valTIMESTEP,ai));
%             fprintf('\tCDi = %0.5f',vecCDI(valTIMESTEP,ai));
        end
        tic
        %% Viscous wrapper
        valCL = vecCL(end,ai);
        valCDI = vecCDI(end,ai);
        
        q_inf = valWEIGHT/(valCL*valAREA);
        valVINF = sqrt(2*q_inf/valDENSITY);
        di = valCDI*valAREA*q_inf;
        
        % Summing freestream and induced forces of each DVE
        % This will NOT work with rotors, it does not take into
        % account freestream! UINF^2*AREA should be the denominator
        vecDVECN = ((vecDVENFREE + vecDVENIND).*2)./(vecDVEAREA);
        vecDVECL = ((vecDVELFREE + vecDVELIND).*2)./(vecDVEAREA);
        vecDVECY = ((vecDVESFREE + vecDVESIND).*2)./(vecDVEAREA);
        
        [ledves, ~, ~] = find(vecDVELE > 0);
        lepanels = vecDVEPANEL(ledves);
        
        vecCNDIST = [];
        vecCLDIST = [];
        vecCYDIST = [];
        matXYZDIST = [];
        matLEDVEDIST = [];
        vecREDIST = [];
        vecAREADIST = [];
        
        for i = 1:max(vecDVEWING)
            
            %% Getting the CL, CY, CN distribution
            idxdve = ledves(vecDVEWING(ledves) == i);
            idxpanel = lepanels(vecDVEWING(ledves) == i);
            
            m = vecM(idxpanel);
            if any(m - m(1))
                disp('Problem with wing chordwise elements.');
                break
            end
            
            m = m(1);
            len = length(vecCLDIST); % start point for this panel in the vectors
            
            % Matrix of how much we need to add to an index to get the next chordwise element
            % It is done this way because n can be different for each panel. Unlike in the wake,
            % we can't just add a constant value to get to the same spanwise location in the next
            % row of elements
            tempm = repmat(vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
            
            rows = repmat(idxdve,1,m) + tempm;
            
            vecCNDIST = [vecCNDIST; sum(vecDVECN(rows),2)];
            vecCLDIST = [vecCLDIST; sum(vecDVECL(rows),2)];
            vecCYDIST = [vecCYDIST; sum(vecDVECY(rows),2)];
            
            % The average coordinates for this row of elements
            matXYZDIST = [matXYZDIST; mean(permute(reshape(matCENTER(rows,:)',3,[],m),[2 1 3]),3)];
            
            % The leading edge DVE for the distribution
            matLEDVEDIST = [matLEDVEDIST; idxdve];
            
            %% Wing/horizontal stabilizer lift and drag
            
            vecREDIST = [vecREDIST; valVINF.*2.*sum(vecDVEHVCRD(rows),2)./valKINV];
            vecAREADIST = [vecAREADIST; sum(vecDVEAREA(rows),2)];
            
            for j = 1:length(idxpanel)
                pan = idxpanel(j);
                airfoil = dlmread(strcat('airfoils/airfoil',num2str(vecAIRFOIL(pan)),'.dat'),'\t', 1, 0);
                
                HiRe = airfoil(end,4);
                LoRe = airfoil(1,4);
                
                cl = vecCLDIST(len + j);
                
                if vecREDIST(len + j) > HiRe
                    if flagVERBOSE == 1
                        disp('Re higher than airfoil Re data')
                    end
                    cl_max = airfoil(end,4);
                elseif vecREDIST(len + j) < LoRe
                    if flagVERBOSE == 1
                        disp('Re lower than airfoil Re data');
                    end
                    Re2 = airfoil(1,4);
                    temp_var = airfoil(airfoil(:,4) == Re2, 2);
                    cl_max = temp_var(end);
                else
                    re1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 4);
                    re1 = re1(end);
                    cl_max1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 2);
                    cl_max1 = cl_max1(end);
                    
                    temp_var = airfoil(airfoil(:,4) > vecREDIST(len + j),4);
                    re2 = temp_var(1);
                    temp_var = airfoil(airfoil(:,4) == (temp_var(1)), 2);
                    cl_max2 = temp_var(end);
                    
                    cl_max = interp1([re1 re2],[cl_max1 cl_max2], vecREDIST(len + j));
                end
                
                % correcting the section cl if we are above cl_max
                if cl > cl_max
                    if flagVERBOSE == 1
                        fprintf('\nStall of Wing %d Section %d, cl = %f Re = %0.0f', i, j, cl, vecREDIST(len + j))
                    end
                    vecCLDIST(len + j) = 0.825*cl_max; % setting the stalled 2d cl
                end
                
                F = scatteredInterpolant(airfoil(:,4), airfoil(:,2), airfoil(:,3),'nearest');
                vecCDPDIST(len + j, 1) = F(vecREDIST(len + j), cl);
                
                % Octave:
                % vecCDPDIST(len + j, 1) = griddata(airfoil(:,4), airfoil(:,2), airfoil(:,3), vecREDIST(len + j), cl, 'nearest');
                
            end
        end
        
        dprof = sum(vecCDPDIST.*q_inf.*vecAREADIST);
        
        % This function does not account for symmetry well, it is all or nothing with symmetry,
        % but it really should be wing-by-wing
        if any(vecSYM) == 1
            dprof = 2.*dprof;
        end
        
        %% Vertical tail drag
        
        dvt = 0;
        for ii = 1:valVSPANELS
            Re = valVINF*matVSGEOM(ii,2)/valKINV;
            
            % Load airfoil data
            airfoil = dlmread(strcat('airfoils/airfoil',num2str(matVSGEOM(ii,4)),'.dat'),'\t', 1, 0);
            
            % determining the drag coefficient corresponding to lift
            % coefficient of 0
            
            % MATLAB:
            F = scatteredInterpolant(airfoil(:,4), airfoil(:,2), airfoil(:,3),'nearest');
            cdvt = F(Re, 0);
            % Octave:
            % cdvt = griddata(Temp.Airfoil(:,4), Temp.Airfoil(:,2), Temp.Airfoil(:,3), Re, 0, 'nearest');
            
            dvt = dvt + cdvt*matVSGEOM(ii,3);
        end
        
        dvt = dvt*q_inf;
        
        %% Fuselage drag
        
        dfuselage = 0;
        
        tempSS = valVINF*valFPWIDTH/valKINV;
        
        for ii = 1:valFPANELS
            Re_fus = (ii-0.5)*tempSS;
            if ii < valFTURB
                cdf = 0.664/sqrt(Re_fus); % Laminar
            else
                cdf = 0.0576/(Re_fus^0.2); % Turbulent
            end
            
            dfuselage = dfuselage + cdf*matFGEOM(ii,2)*pi*valFPWIDTH;
        end
        
        dfuselage = dfuselage*q_inf;
        
        %% Total Drag
        
        valD = di + dprof + dvt + dfuselage;
        
        dint = valD*(valINTERF/100);
        
        valD = valD + dint;
        
        
        
        
        
        toc
    end
end

fprintf('\n');

%% Plotting

% if flagPLOT == 1
%     [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
%     [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
%     [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');

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

% end

%% Viscous wrapper

% whos