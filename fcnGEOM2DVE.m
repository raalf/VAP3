function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, OUTP, SURF)

% tranlsate INPU.matGEOM to vehicle origin
INPU.matGEOM(:,1:3,:) = INPU.matGEOM(:,1:3,:)+permute(reshape(INPU.matVEHORIG(INPU.matGEOM(:,6,:),:)',3,2,[]),[2,1,3]);

vecPANELSURFACE = sort([INPU.vecPANELWING INPU.vecPANELROTOR + max(INPU.vecPANELWING).*uint16((INPU.vecPANELROTOR > 0))],2,'descend');
vecPANELSURFACE = vecPANELSURFACE(:,1);

SURF.valWINGS = max(vecPANELSURFACE);

SURF.matCENTER = [];
SURF.vecDVEHVSPN = [];
SURF.vecDVEHVCRD = [];
SURF.vecDVELESWP = [];
SURF.vecDVEMCSWP = [];
SURF.vecDVETESWP = [];
SURF.vecDVEROLL = [];
SURF.vecDVEPITCH = [];
SURF.vecDVEYAW = [];
SURF.vecDVEAREA = [];
SURF.matDVENORM = [];
SURF.matVLST = [];
SURF.matNTVLST = [];
SURF.matNPVLST = [];
SURF.matDVE = [];
SURF.valNELE = 0;
SURF.matADJE = [];
SURF.vecDVESYM = [];
SURF.vecDVETIP = [];
SURF.vecDVESURFACE = [];
SURF.vecDVELE = [];
SURF.vecDVETE = [];
SURF.vecDVEPANEL = [];
SURF.matPANELTE = [];
SURF.vecDVETRI = [];
SURF.vecWINGTYPE = [];

INPU.vecPANELROTOR = uint8(INPU.vecPANELROTOR);

WAKE.valWSIZETRI = 0;
WAKE.valWSIZE = 0;

for i = unique(vecPANELSURFACE,'stable')'
    
    panels = length(nonzeros(vecPANELSURFACE == i));
   

        [tSURF.matCENTER, tSURF.vecDVEHVSPN, tSURF.vecDVEHVCRD, tSURF.vecDVELESWP, tSURF.vecDVEMCSWP, tSURF.vecDVETESWP, ...
            tSURF.vecDVEROLL, tSURF.vecDVEPITCH, tSURF.vecDVEYAW, tSURF.vecDVEAREA, tSURF.matDVENORM, ...
            tSURF.matVLST, tSURF.matNTVLST, tSURF.matNPVLST, tSURF.matDVE, tSURF.valNELE, tSURF.matADJE, ...
            tSURF.vecDVESYM, tSURF.vecDVETIP, tSURF.vecDVESURFACE, tSURF.vecDVELE, tSURF.vecDVETE, tSURF.vecDVEPANEL, tSURF.matPANELTE] = ...
            fcnGENERATEDVES(panels, INPU.matGEOM(:,:,(vecPANELSURFACE == i)), INPU.vecSYM(vecPANELSURFACE == i), INPU.vecN(vecPANELSURFACE == i), INPU.vecM(vecPANELSURFACE == i));
        
        SURF.vecDVETRI = [SURF.vecDVETRI; zeros(tSURF.valNELE,1)];
        
    SURF.vecWINGTYPE = [SURF.vecWINGTYPE; repmat(INPU.vecWINGTYPE(i),size(tSURF.matCENTER,1),1)];
        
    SURF.valNELE = SURF.valNELE + tSURF.valNELE;
    SURF.matCENTER = [SURF.matCENTER; tSURF.matCENTER];
    SURF.vecDVEHVSPN = [SURF.vecDVEHVSPN; tSURF.vecDVEHVSPN];
    SURF.vecDVEHVCRD = [SURF.vecDVEHVCRD; tSURF.vecDVEHVCRD];
    SURF.vecDVELESWP = [SURF.vecDVELESWP; tSURF.vecDVELESWP];
    SURF.vecDVEMCSWP = [SURF.vecDVEMCSWP; tSURF.vecDVEMCSWP];
    SURF.vecDVETESWP = [SURF.vecDVETESWP; tSURF.vecDVETESWP];
    SURF.vecDVEROLL = [SURF.vecDVEROLL; tSURF.vecDVEROLL];
    SURF.vecDVEPITCH = [SURF.vecDVEPITCH; tSURF.vecDVEPITCH];
    SURF.vecDVEYAW = [SURF.vecDVEYAW; tSURF.vecDVEYAW];
    SURF.vecDVEAREA = [SURF.vecDVEAREA; tSURF.vecDVEAREA];
    SURF.matDVENORM = [SURF.matDVENORM; tSURF.matDVENORM];
    SURF.vecDVESYM = [SURF.vecDVESYM; tSURF.vecDVESYM];
    SURF.vecDVETIP = [SURF.vecDVETIP; tSURF.vecDVETIP];
    SURF.vecDVELE = [SURF.vecDVELE; tSURF.vecDVELE];
    SURF.vecDVETE = [SURF.vecDVETE; tSURF.vecDVETE];
    SURF.matPANELTE = [SURF.matPANELTE; tSURF.matPANELTE];
    
    if i == 1; surfaceoffset = 0; paneloffset = 0;
    else; surfaceoffset = max(SURF.vecDVESURFACE); paneloffset = max(SURF.vecDVEPANEL);
    end
    
    SURF.vecDVESURFACE = [SURF.vecDVESURFACE; uint8(uint8(tSURF.vecDVESURFACE) + surfaceoffset)];
    SURF.vecDVEPANEL = [SURF.vecDVEPANEL; uint16(uint16(tSURF.vecDVEPANEL) + paneloffset)];
    
    vlstoffset = size(SURF.matVLST,1);
    dveoffset = size(SURF.matDVE,1);
    SURF.matDVE = [SURF.matDVE; tSURF.matDVE + vlstoffset];
    SURF.matVLST = [SURF.matVLST; tSURF.matVLST];
    SURF.matNTVLST = [SURF.matNTVLST; tSURF.matNTVLST];
    SURF.matNPVLST = [SURF.matNPVLST; tSURF.matNPVLST];
  
    tSURF.matADJE = [tSURF.matADJE(:,1) + dveoffset tSURF.matADJE(:,2) tSURF.matADJE(:,3) + dveoffset tSURF.matADJE(:,4)];
    SURF.matADJE = [SURF.matADJE; tSURF.matADJE];
end

% % Identifying which DVEs belong to which vehicle, as well as which type of lifting surface they belong to (wing or rotor)
SURF.vecDVEVEHICLE = VEHI.vecSURFACEVEHICLE(SURF.vecDVESURFACE);
SURF.vecDVEWING = SURF.vecDVESURFACE;

SURF.vecDVEROTOR = INPU.vecPANELROTOR(SURF.vecDVEPANEL); % Alton-Y
SURF.vecDVEROTORBLADE = double(any(SURF.vecDVEROTOR,2)); % Current rotor DVEs are for Blade 1 (they are duplicated to Blade 2, 3, etc etc below)
idx_rotor = SURF.vecDVEROTOR>0; % Alton-Y
SURF.vecDVEWING(idx_rotor) = 0;

MISC.matSURFACETYPE = uint8(zeros(size(unique(SURF.vecDVESURFACE),1),2));
MISC.matSURFACETYPE(nonzeros(unique(SURF.vecDVEWING)),1) = nonzeros(unique(SURF.vecDVEWING));
MISC.matSURFACETYPE(nonzeros(unique(SURF.vecDVESURFACE(idx_rotor))),2) = nonzeros(unique(SURF.vecDVEROTOR));

% Identifying which ROTOR belongs to which vehicle.
VEHI.vecROTORVEH = VEHI.vecSURFACEVEHICLE(MISC.matSURFACETYPE(:,2)~=0);

% Duplicate Blades in a Rotor
[INPU.vecPANELROTOR, INPU.vecN, INPU.vecM, SURF.matVLST, SURF.matCENTER, SURF.matDVE, SURF.matADJE, SURF.vecDVEVEHICLE, ...
    SURF.vecDVEWING, SURF.vecDVEROTOR, MISC.matSURFACETYPE, SURF.vecDVESURFACE, SURF.vecDVEPANEL, ...
    SURF.vecDVETIP, SURF.vecDVELE, SURF.vecDVETE, SURF.vecDVEROTORBLADE, SURF.vecDVESYM, ...
    SURF.valNELE, SURF.matNTVLST, VISC.cellAIRFOIL, OUTP.ROTOR] = fcnDUPBLADE( VEHI.vecROTORVEH, SURF.vecDVEROTOR, ...
    SURF.matVLST, SURF.matCENTER, SURF.matDVE, SURF.matADJE, INPU.vecROTORBLADES, ...
    SURF.valNELE, INPU.matROTORHUB, INPU.matVEHORIG, SURF.vecDVEVEHICLE, SURF.vecDVEWING, ...
    MISC.matSURFACETYPE, SURF.vecDVESURFACE, SURF.vecDVEPANEL, SURF.vecDVETIP, SURF.vecDVELE, ...
    SURF.vecDVETE, SURF.vecDVEROTORBLADE, SURF.vecDVESYM, INPU.matROTORAXIS, SURF.matNTVLST, ...
    INPU.vecM, INPU.vecN, INPU.vecPANELROTOR, VISC.cellAIRFOIL);

VEHI.vecBFRAME = [-1;0;0]; % Unit vector for vehicle body frame direction

[ VEHI.matGLOBUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, COND.vecVEHALPHA, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
COND.vecVEHPITCH = rad2deg(VEHI.matVEHROT(2));
[SURF.matVLST, SURF.matCENTER, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNPVLST, VEHI.vecBFRAME] = fcnROTVEHICLE( SURF.matDVE, SURF.matVLST, SURF.matCENTER, INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH, SURF.matNTVLST, VEHI.vecBFRAME);

[ SURF.matUINF ] = fcnINITUINF( SURF.matCENTER, VEHI.matGLOBUVW, VEHI.matVEHROT, SURF.vecDVEVEHICLE, ...
    SURF.vecDVEROTOR, VEHI.vecROTORVEH, INPU.matVEHORIG, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, COND.vecROTORRPM );

% update DVE params after vehicle rotation
[ SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW,...
    SURF.vecDVELESWP, SURF.vecDVEMCSWP, SURF.vecDVETESWP, SURF.vecDVEAREA, SURF.matDVENORM, ~, ~, SURF.matCENTER ] ...
    = fcnVLST2DVEPARAM(SURF.matDVE, SURF.matVLST);

WAKE.valWSIZE = length(nonzeros(SURF.vecDVETE));

% Compute torque arm length for rotor power calculations
SURF.vecQARM = zeros(SURF.valNELE,3);
if max(SURF.vecDVEROTOR)>0
    SURF.vecQARM(SURF.vecDVEROTOR>0,:) = SURF.matCENTER(SURF.vecDVEROTOR>0,:) - INPU.matROTORHUB(SURF.vecDVEROTOR(SURF.vecDVEROTOR>0),:);    
    SURF.vecQARM = sqrt(sum(SURF.vecQARM.^2,2));
end

% Compute pitch arm for tail if it exists
if any(SURF.vecWINGTYPE == 2)
    SURF.vecPITCHARM = [];

    for j = unique(SURF.vecWINGTYPE,'stable')'

        [ledves, ~, ~] = find(SURF.vecDVELE > 0);
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

        vecCRDDIST(isCurWing,1) = sum(2*SURF.vecDVEHVCRD(rows),2);

        % Find LE mid-pt location in xyz of each DVE
        temp = fcnGLOBSTAR(SURF.matCENTER(ledves(isCurWing),:), SURF.vecDVEROLL(ledves(isCurWing)), SURF.vecDVEPITCH(ledves(isCurWing)), SURF.vecDVEYAW(ledves(isCurWing)));
        SURF.matDVEQTRCRD = fcnSTARGLOB([temp(:,1)-SURF.vecDVEHVCRD(ledves(isCurWing)) + 0.25*vecCRDDIST(isCurWing,1),temp(:,2),temp(:,3)], SURF.vecDVEROLL(ledves(isCurWing)), SURF.vecDVEPITCH(ledves(isCurWing)), SURF.vecDVEYAW(ledves(isCurWing)));
        
        SURF.vecPITCHARM = [SURF.vecPITCHARM; (INPU.vecVEHCG(1,1) - SURF.matDVEQTRCRD(:,1)).*cos(pi*COND.vecVEHALPHA(1)/180)...
            + (INPU.vecVEHCG(1,3) - SURF.matDVEQTRCRD(:,3)).*sin(pi*COND.vecVEHALPHA(1)/180)];
        
        SURF.matLEMIDPT = (SURF.matVLST(SURF.matDVE(rows,1),:) + SURF.matVLST(SURF.matDVE(rows,2),:))./2;

        % Compute pitch moment arm from CG to DVE LE midpoint
%         SURF.vecPITCHARM = [SURF.vecPITCHARM; (INPU.vecVEHCG(1,1) - SURF.matLEMIDPT(:,1)).*cos(pi*COND.vecVEHALPHA/180)...
%                 + (INPU.vecVEHCG(1,3) - SURF.matLEMIDPT(:,3)).*sin(pi*COND.vecVEHALPHA/180)];
    end
end

end
