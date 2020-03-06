function [SURF, OUTP] = fcnPITCHMOMENT(FLAG, SURF, OUTP, INPU, COND, VEHI)
% #tbt to this function. It's baaaccccccckkkkkkkkkk
% Computing the total pitching moment for the vehicle

OUTP.vecVEHPITCHMOM = [];
SURF.matDVEQTRCRD = [];
matLIFTDIR = [];
matDRAGDIR = [];
pitcharm = [];
deltaM_lift = [];
deltaM_drag = [];
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
    if any(SURF.vecWINGTYPE == 2)
        SURF.vecPITCHARM = [];

        vecCRDDIST(isCurWing,1) = sum(2*SURF.vecDVEHVCRD(rows),2);

        % Find LE mid-pt location in xyz of each DVE
        temp = fcnGLOBSTAR(SURF.matCENTER(ledves(isCurWing),:), SURF.vecDVEROLL(ledves(isCurWing)), SURF.vecDVEPITCH(ledves(isCurWing)), SURF.vecDVEYAW(ledves(isCurWing)));
        matDVEQTRCRD = fcnSTARGLOB([temp(:,1)-SURF.vecDVEHVCRD(ledves(isCurWing)) + 0.25*vecCRDDIST(isCurWing,1),temp(:,2),temp(:,3)], SURF.vecDVEROLL(ledves(isCurWing)), SURF.vecDVEPITCH(ledves(isCurWing)), SURF.vecDVEYAW(ledves(isCurWing)));
        
        lemidpt = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2;
        
        SURF.matDVEQTRCRD = [SURF.matDVEQTRCRD; matDVEQTRCRD];
        matLIFTDIR = [matLIFTDIR; SURF.matLIFTDIR(ledves(isCurWing),:)];
        matDRAGDIR = [matDRAGDIR; SURF.matDRAGDIR(isCurWing,:)];
        SURF.vecPITCHARM = SURF.matDVEQTRCRD - INPU.vecVEHCG;
        
        pitcharm = [pitcharm; lemidpt - INPU.vecVEHCG];
        
%         SURF.vecPITCHARM = fcnGLOBSTAR(SURF.vecPITCHARM, deg2rad(COND.vecVEHROLL*ones(size(SURF.vecPITCHARM,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(SURF.vecPITCHARM,1),1)), deg2rad(COND.vecVEHBETA*ones(size(SURF.vecPITCHARM,1),1)));
        
    end

        % Compute pitch moment from lifting lines
        if FLAG.VISCOUS == 1
            [dp_wing, dp_tail] = SURF.vecDPDIST.DPDIST; 
            dpdist = [dp_wing; dp_tail]; % Visocus drag force distribution
        else
            dpdist = zeros(size(SURF.vecDVEDIND(tedves,:),1),1);
            OUTP.vecCMDIST = zeros(length(isCurWing),1);
        end
        
        if isempty(VEHI.vecPROPLOC) == 0 && j > 1
            thrust = sum(sum(SURF.vecDVEDIND(tedves,:).*matDRAGDIR + dpdist.*matDRAGDIR,1),2).*VEHI.vecPROPDIR;
            deltaM_prop = cross(VEHI.vecPROPLOC - INPU.vecVEHCG, thrust);
        else
            deltaM_prop = [0 0 0];
        end
        
        lift = (SURF.vecDVELFREE(SURF.vecWINGTYPE == j) + SURF.vecDVELIND(SURF.vecWINGTYPE == j)).*SURF.matLIFTDIR(SURF.vecWINGTYPE == j,:);
        drag = (SURF.vecDVEDIND(SURF.vecWINGTYPE(tedves) == j)).*SURF.matDRAGDIR(isCurWing,:);
        
%         deltaM_lift = [deltaM_lift; cross(pitcharm((SURF.vecWINGTYPE == j),:),lift)];
%         deltaM_drag = [deltaM_drag; cross(pitcharm(tedves(isCurWing),:),drag)];
        
        tempforce = (sum(SURF.vecDVELFREE(rows),2) + sum(SURF.vecDVELIND(rows),2)).*matLIFTDIR(isCurWing,:) + SURF.vecDVEDIND(tedves(isCurWing),:).*matDRAGDIR(isCurWing,:) + dpdist(isCurWing).*matDRAGDIR(isCurWing,:);
        deltaM = cross(SURF.vecPITCHARM(isCurWing,:),tempforce);
        
%         deltaM = sum(deltaM_lift(:,2),1) + sum(deltaM_drag(:,2),1);
        
        % Total vehicle pitch moment (Moment/density)
        
        if j == 2
            OUTP.vecCMDIST(isCurWing) = OUTP.vecCMDIST(isCurWing)*0;
        end
        
        if any(isnan(OUTP.vecCMDIST)) == 1
            OUTP.vecCMDIST(isnan(OUTP.vecCMDIST)) = 0;
        end
        
        % Double vehicle pitch moment for symmetry
        if any(INPU.vecSYM == 1) == 1
            OUTP.vecVEHPITCHMOM(i,j) = 2*(sum(deltaM(:,2),1) + sum(OUTP.vecCMDIST(isCurWing).*(0.5*COND.vecVEHVINF*COND.vecVEHVINF*vecAREADIST(isCurWing)'*INPU.vecCMAC(i)),1));
%             OUTP.vecVEHPITCHMOM(i,j) = 2*(deltaM + sum(OUTP.vecCMDIST(isCurWing).*(0.5*COND.vecVEHVINF*COND.vecVEHVINF*vecAREADIST(isCurWing)'*INPU.vecCMAC(i)),1));
        else
            OUTP.vecVEHPITCHMOM(i,j) = sum(deltaM,1);
        end
    
    end
    
    OUTP.vecVEHPITCHMOM = sum(OUTP.vecVEHPITCHMOM(i,:),2) + 2*deltaM_prop(:,2);
    
    % Non-dimensionalize
    OUTP.vecVEHCM(i) = OUTP.vecVEHPITCHMOM/(0.5*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA(i)*INPU.vecCMAC(i));

end