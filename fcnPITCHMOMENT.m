function [SURF, OUTP] = fcnPITCHMOMENT(FLAG, SURF, OUTP, INPU, COND, VEHI)
% #tbt to this function. It's baaaccccccckkkkkkkkkk
% Computing the total pitching moment for the vehicle

OUTP.vecVEHPITCHMOM = [];
SURF.matDVEQTRCRD = [];
matLIFTDIR = [];
matDRAGDIR = [];
pitcharm = [];
deltaM_norm = [];
deltaM_drag = [];

%% Aerodynamic center method
% Find perpendicular distance from DVE control point to CG for each vehicle
% for i = 1:INPU.valVEHICLES
%     
%     for j = unique(SURF.vecWINGTYPE,'stable')'
%         
%         if j == 3
%             break;
%         end
%         
%         [ledves, ~, ~] = find(SURF.vecDVELE > 0);
%         [tedves, ~, ~] = find(SURF.vecDVETE > 0);
%         lepanels = SURF.vecDVEPANEL(ledves);
%         
%         isCurWing = SURF.vecWINGTYPE(ledves) == j;
%         
%         idxdve = uint16(ledves(isCurWing));
%         idxpanel = lepanels(isCurWing);
% 
% 
%         m = INPU.vecM(idxpanel);
%         if any(m - m(1))
%             disp('Problem with wing chordwise elements.');
%             break
%         end
% 
%         m = m(1);
% 
%         % Matrix of how much we need to add to an index to get the next chordwise element
%         % It is done this way because n can be different for each panel. Unlike in the wake,
%         % we can't just add a constant value to get to the same spanwise location in the next
%         % row of elements
%         tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);
% 
%         rows = repmat(idxdve,1,m) + uint16(tempm);
%         
%         vecAREADIST(isCurWing) = sum(SURF.vecDVEAREA(rows),2);
%         
%         % Compute pitch arm for tail if it exists
%     if any(SURF.vecWINGTYPE == 2)
%         SURF.vecPITCHARM = [];
% 
%         vecCRDDIST(isCurWing,1) = sum(2*SURF.vecDVEHVCRD(rows),2);
%         
%         lemidpt = (SURF.matVLST(SURF.matDVE(rows(:,1),1),:) + SURF.matVLST(SURF.matDVE(rows(:,1),2),:))./2;
%         lemidpt = fcnGLOBSTAR(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE(rows(:,1)) == j),SURF.vecDVEPITCH(SURF.vecWINGTYPE(rows(:,1)) == j),SURF.vecDVEYAW(SURF.vecWINGTYPE(rows(:,1)) == j));
%         lemidpt(:,1) = lemidpt(:,1) + 0.25*vecCRDDIST(isCurWing);
%         vecQTRCRD = fcnSTARGLOB(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE(rows(:,1)) == j),SURF.vecDVEPITCH(SURF.vecWINGTYPE(rows(:,1)) == j),SURF.vecDVEYAW(SURF.vecWINGTYPE(rows(:,1)) == j));
% 
%         matLIFTDIR = [matLIFTDIR; SURF.matLIFTDIR(ledves(isCurWing),:)];
%         matDRAGDIR = [matDRAGDIR; SURF.matDRAGDIR(isCurWing,:)];
%         SURF.vecPITCHARM = vecQTRCRD - INPU.vecVEHCG; % Pitch arm to CG
% %         SURF.vecPITCHARM = SURF.matDVEQTRCRD - SURF.matEALST(1,:); % Pitch arm to wing elastic axis
%         
%         % Convert pitch arm to body frame
%         SURF.vecPITCHARM = fcnGLOBSTAR(SURF.vecPITCHARM, deg2rad(COND.vecVEHROLL*ones(size(SURF.vecPITCHARM,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(SURF.vecPITCHARM,1),1)), deg2rad(COND.vecVEHBETA*ones(size(SURF.vecPITCHARM,1),1)));
%         
%     end
% 
%         % Compute pitch moment from lifting lines
%         if FLAG.VISCOUS == 1
%             [dp_wing, dp_tail] = SURF.vecDPDIST.DPDIST; 
%             dpdist = [dp_wing; dp_tail]; % Visocus drag force distribution
%         else
%             dpdist = zeros(size(SURF.vecDVEDIND(tedves,:),1),1);
%             OUTP.vecCMDIST = zeros(length(isCurWing),1);
%         end
%         
%         % Compute total force vector and pitching moment due to it. Moment
%         % is per density
%         tot_force = (sum(SURF.vecDVENFREE(rows),2) + sum(SURF.vecDVENIND(rows),2)).*SURF.matDVENORM(ledves(isCurWing),:) + SURF.vecDVEDIND(tedves(isCurWing),:).*matDRAGDIR(isCurWing,:);
%         deltaM = cross(SURF.vecPITCHARM,tot_force);
%         deltaM = sum(deltaM(:,2),1);
%         
%         if j == 2
%             OUTP.vecCMDIST(isCurWing) = OUTP.vecCMDIST(isCurWing)*0;
%         end
%         
%         if any(isnan(OUTP.vecCMDIST)) == 1
%             OUTP.vecCMDIST(isnan(OUTP.vecCMDIST)) = 0;
%         end
%         
%         % Double vehicle pitch moment for symmetry
%         if any(INPU.vecSYM == 1) == 1
%             OUTP.vecVEHPITCHMOM(i,j) = 2*(deltaM + sum(OUTP.vecCMDIST(isCurWing).*(0.5*COND.vecVEHVINF*COND.vecVEHVINF*vecAREADIST(isCurWing)'*INPU.vecCMAC(i)),1));
%         else
%             OUTP.vecVEHPITCHMOM(i,j) = sum(deltaM,1);
%         end
%     
%     end
%     
%     OUTP.vecVEHPITCHMOM = sum(OUTP.vecVEHPITCHMOM(i,:),2).*COND.valDENSITY;
%     
%     % Non-dimensionalize
%     OUTP.vecVEHCM(i) = OUTP.vecVEHPITCHMOM/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA(i)*INPU.vecCMAC(i));
% 
% end

%% Lifting line method
% Find perpendicular distance from DVE control point to CG for each vehicle
for i = 1:INPU.valVEHICLES
    
    for j = unique(SURF.vecWINGTYPE,'stable')'
        
        if j == 3
            break;
        end
        
        [ledves, ~, ~] = find(SURF.vecDVELE > 0);
        [tedves, ~, ~] = find(SURF.vecDVETE > 0);
        lepanels = SURF.vecDVEPANEL(ledves);
        
        isCurWing = SURF.vecWINGTYPE(ledves) == j;
        
        idxdve = uint16(ledves(isCurWing));
        idxpanel = lepanels(isCurWing);


        m = INPU.vecM(idxpanel);
        if any(m - m(1))
            disp('Problem with wing chordwise elements.');
            break
        end

        m = m(1);

        % Matrix of how much we need to add to an index to get the next chordwise element
        % It is done this way because n can be different for each panel. Unlike in the wake,
        % we can't just add a constant value to get to the same spanwise location in the next
        % row of elements
        tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);

        rows = repmat(idxdve,1,m) + uint16(tempm);
        
        vecAREADIST(isCurWing) = sum(SURF.vecDVEAREA(rows),2);
        
        % Compute pitch arm for tail if it exists
%     if any(SURF.vecWINGTYPE == 2)
        SURF.vecPITCHARM = [];

        vecCRDDIST(isCurWing,1) = sum(2*SURF.vecDVEHVCRD(rows),2);

        % Find LE mid-pt location in xyz of each DVE
        lemidpt = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2;
        lemidpt = fcnGLOBSTAR(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == j),SURF.vecDVEPITCH(SURF.vecWINGTYPE == j),SURF.vecDVEYAW(SURF.vecWINGTYPE == j));
        lemidpt(:,1) = lemidpt(:,1) + SURF.vecDVEHVCRD(SURF.vecWINGTYPE == j)./2;
        lemidpt = fcnSTARGLOB(lemidpt,SURF.vecDVEROLL(SURF.vecWINGTYPE == j),SURF.vecDVEPITCH(SURF.vecWINGTYPE == j),SURF.vecDVEYAW(SURF.vecWINGTYPE == j));
%         matLIFTDIR = [matLIFTDIR; SURF.matLIFTDIR(ledves(isCurWing),:)];
        matDRAGDIR = [matDRAGDIR; SURF.matDRAGDIR(isCurWing,:)];
        
        if j == 1
%             temp_pitcharm = lemidpt - SURF.matEALST(rows,:);
        else
%             temp_pitcharm = lemidpt - SURF.matEALST(1,:);
        end
        temp_pitcharm = lemidpt - INPU.vecVEHCG;
%         temp_pitcharm = fcnGLOBSTAR(temp_pitcharm, deg2rad(COND.vecVEHROLL*ones(size(temp_pitcharm,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(temp_pitcharm,1),1)), deg2rad(COND.vecVEHBETA*ones(size(temp_pitcharm,1),1)));
        pitcharm = [pitcharm; temp_pitcharm];
                
%     end

        % Compute pitch moment from lifting lines
        if FLAG.VISCOUS == 1
            [dp_wing, dp_tail] = SURF.vecDPDIST.DPDIST; 
            dpdist = [dp_wing; dp_tail]; % Visocus drag force distribution
        else
            dpdist = zeros(size(SURF.vecDVEDIND(tedves,:),1),1);
            OUTP.vecCMDIST = zeros(length(isCurWing),1);
        end
        
        % Moment due to normal force and drag
%         norm = (SURF.vecDVENFREE(SURF.vecWINGTYPE == j) + SURF.vecDVENIND(SURF.vecWINGTYPE == j)).*SURF.matNORMDIR(SURF.vecWINGTYPE == j,:);
%         drag = (SURF.vecDVEDIND(tedves(isCurWing))).*SURF.matDRAGDIR(isCurWing,:);

%         lift = dot(SURF.matDVEIFORCE,repmat(VEHI.ldir,size(SURF.matDVEIFORCE,1),1),2).*VEHI.ldir;
%         drag = (SURF.vecDVEDIND(tedves(isCurWing))).*SURF.matDRAGDIR(isCurWing,:);
        lift = dot(SURF.matDVEIFORCE,VEHI.ldir,2).*VEHI.ldir;
        drag = (SURF.vecDVEDIND(tedves(isCurWing))).*SURF.matDRAGDIR(isCurWing,:);
        
%         deltaM_norm = cross(pitcharm((SURF.vecWINGTYPE == j),:),norm);
        deltaM_norm = cross(pitcharm((SURF.vecWINGTYPE == j),:),lift(SURF.vecWINGTYPE == j,:));
        deltaM_drag = cross(pitcharm(tedves(isCurWing),:),drag);
        
        % Total vehicle pitch moment (Moment/density)
        deltaM = sum(deltaM_norm(:,2),1) + sum(deltaM_drag(:,2),1);
        
         % Moment due to point masses. Only used if not computing moment
         % about CG
%         deltaM_mass = cross(VEHI.vecFUSECG-SURF.matEALST(1,:),[0,0,-VEHI.vecFUSEMASS*9.81]) + sum(cross(VEHI.vecWINGCG(2:end,:)-SURF.matEALST(1,:),[zeros(size(VEHI.vecWINGMASS(2:end),1),2),-VEHI.vecWINGMASS(2:end)*9.81]),1); % Pitching moment due to masses of wing, fuse, tail, etc.      
        
        if j == 2
            OUTP.vecCMDIST(isCurWing) = OUTP.vecCMDIST(isCurWing)*0;
        end
        
        if any(isnan(OUTP.vecCMDIST)) == 1
            OUTP.vecCMDIST(isnan(OUTP.vecCMDIST)) = 0;
        end
        
        % Double vehicle pitch moment for symmetry
        if any(INPU.vecSYM == 1) == 1
            OUTP.vecVEHPITCHMOM(i,j) = 2*(deltaM + sum(OUTP.vecCMDIST(isCurWing).*(0.5*COND.vecVEHVINF*COND.vecVEHVINF*vecAREADIST(isCurWing)'*INPU.vecCMAC(i)),1));
        else
            OUTP.vecVEHPITCHMOM(i,j) = sum(deltaM,1);
        end
    
    end
    
    OUTP.vecVEHPITCHMOM = sum(OUTP.vecVEHPITCHMOM(i,:),2).*COND.valDENSITY;
    
    % Non-dimensionalize
    OUTP.vecVEHCM(i) = OUTP.vecVEHPITCHMOM/(0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA(i)*INPU.vecCMAC(i));

end