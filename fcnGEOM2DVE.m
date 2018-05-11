function [INPU, COND, MISC, VISC, WAKE, VEHI, SURF, OUTP] = fcnGEOM2DVE(INPU, COND, VISC, VEHI, WAKE, FLAG, OUTP)

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
SURF.matDVE = uint16([]);
SURF.valNELE = 0;
SURF.matADJE = [];
SURF.vecDVESYM = uint8([]);
SURF.vecDVETIP = uint8([]);
SURF.vecDVESURFACE = uint8([]);
SURF.vecDVELE = uint8([]);
SURF.vecDVETE = uint8([]);
SURF.vecDVEPANEL = uint16([]);
SURF.matPANELTE = uint16([]);
SURF.vecDVETRI = [];

INPU.vecPANELROTOR = uint8(INPU.vecPANELROTOR);

WAKE.valWSIZETRI = 0;
WAKE.valWSIZE = 0;

for i = unique(vecPANELSURFACE,'stable')'
    
    panels = length(nonzeros(vecPANELSURFACE == i));
   
%     if isnan(COND.vecWINGTRI(i))
        [tSURF.matCENTER, tSURF.vecDVEHVSPN, tSURF.vecDVEHVCRD, tSURF.vecDVELESWP, tSURF.vecDVEMCSWP, tSURF.vecDVETESWP, ...
            tSURF.vecDVEROLL, tSURF.vecDVEPITCH, tSURF.vecDVEYAW, tSURF.vecDVEAREA, tSURF.matDVENORM, ...
            tSURF.matVLST, tSURF.matNTVLST, tSURF.matDVE, tSURF.valNELE, tSURF.matADJE, ...
            tSURF.vecDVESYM, tSURF.vecDVETIP, tSURF.vecDVESURFACE, tSURF.vecDVELE, tSURF.vecDVETE, tSURF.vecDVEPANEL, tSURF.matPANELTE] = ...
            fcnGENERATEDVES(panels, INPU.matGEOM(:,:,(vecPANELSURFACE == i)), INPU.vecSYM(vecPANELSURFACE == i), INPU.vecN(vecPANELSURFACE == i), INPU.vecM(vecPANELSURFACE == i));
        
        SURF.vecDVETRI = [SURF.vecDVETRI; zeros(tSURF.valNELE,1)];
        
%     else
%         
%         [tSURF.matCENTER, tSURF.vecDVEHVSPN, tSURF.vecDVEHVCRD, tSURF.vecDVELESWP, tSURF.vecDVEMCSWP, tSURF.vecDVETESWP, ...
%             tSURF.vecDVEROLL, tSURF.vecDVEPITCH, tSURF.vecDVEYAW, tSURF.vecDVEAREA, tSURF.matDVENORM, ...
%             tSURF.matVLST, tSURF.matNTVLST, tSURF.matDVE, tSURF.valNELE, tSURF.matADJE, ...
%             tSURF.vecDVESYM, tSURF.vecDVETIP, tSURF.vecDVESURFACE, tSURF.vecDVELE, tSURF.vecDVETE, ...
%             tSURF.vecDVEPANEL, tmatNPVLST0, INPU.vecM(vecPANELSURFACE == i,1), INPU.vecN(vecPANELSURFACE == i), tSURF.matPANELTE] = ...
%             fcnGENERATEDVESTRI(panels, INPU.matGEOM(:,:,(vecPANELSURFACE == i)), INPU.vecSYM(vecPANELSURFACE == i), INPU.vecN(vecPANELSURFACE == i), INPU.vecM(vecPANELSURFACE == i));
%         
%         SURF.vecDVETRI = [SURF.vecDVETRI; ones(tSURF.valNELE,1)];
        
        
%     end
    
%     if isnan(COND.vecWAKETRI(i))

%     else
%         WAKE.valWSIZETRI = WAKE.valWSIZETRI + sum(nonzeros(tSURF.vecDVETE > 0))*2;
%     end
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
  
    tSURF.matADJE = [tSURF.matADJE(:,1) + dveoffset tSURF.matADJE(:,2) tSURF.matADJE(:,3) + dveoffset tSURF.matADJE(:,4)];
    SURF.matADJE = [SURF.matADJE; tSURF.matADJE];
    temp = INPU.vecM(vecPANELSURFACE == i).*INPU.vecN(vecPANELSURFACE == i);
% 	tempM = [tempM; repmat(INPU.vecM(vecPANELSURFACE == i),[temp(1),1])];
%     tempN = [tempN; repmat(INPU.vecN(vecPANELSURFACE == i),[temp(1),1])];
end

% INPU.vecM = tempM;
% INPU.vecN = tempN;

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
    SURF.valNELE, SURF.matNTVLST, VISC.cellAIRFOIL, OUTP.DISTRIBUTIONS.matBLADEXYPLANE] = fcnDUPBLADE( VEHI.vecROTORVEH, SURF.vecDVEROTOR, ...
    SURF.matVLST, SURF.matCENTER, SURF.matDVE, SURF.matADJE, INPU.vecROTORBLADES, ...
    SURF.valNELE, INPU.matROTORHUB, INPU.matVEHORIG, SURF.vecDVEVEHICLE, SURF.vecDVEWING, ...
    MISC.matSURFACETYPE, SURF.vecDVESURFACE, SURF.vecDVEPANEL, SURF.vecDVETIP, SURF.vecDVELE, ...
    SURF.vecDVETE, SURF.vecDVEROTORBLADE, SURF.vecDVESYM, INPU.matROTORAXIS, SURF.matNTVLST, ...
    INPU.vecM, INPU.vecN, INPU.vecPANELROTOR, VISC.cellAIRFOIL);

    if ~isempty(VISC.matFDVE)
        if FLAG.PRINT == 1
            disp('FUSELAGE NOT PROPERLY IMPLEMENTED');
        end

        VISC.matFVLST = VISC.matFVLST + VISC.matFUSEORIG;
    
    end

% flap = 40;
% idx = find(SURF.vecDVETE > 0 & SURF.vecDVEWING > 0);
% for i = idx'
% u = SURF.matVLST(SURF.matDVE(i,2),:) - SURF.matVLST(SURF.matDVE(i,1),:);
% uo = (SURF.matVLST(SURF.matDVE(i,2),:) + SURF.matVLST(SURF.matDVE(i,1),:))./2;
% 
% R = columbia_rotation(u, -flap);
% SURF.matVLST(SURF.matDVE(i,3),:) = (SURF.matVLST(SURF.matDVE(i,3),:) - uo)*R + uo;
% SURF.matVLST(SURF.matDVE(i,4),:) = (SURF.matVLST(SURF.matDVE(i,4),:) - uo)*R + uo;
% end



[ VEHI.matVEHUVW, VEHI.matVEHROT, VEHI.matVEHROTRATE, MISC.matCIRORIG] = fcnINITVEHICLE( COND.vecVEHVINF, INPU.matVEHORIG, COND.vecVEHALPHA, COND.vecVEHBETA, COND.vecVEHFPA, COND.vecVEHROLL, COND.vecVEHTRK, VEHI.vecVEHRADIUS );
[SURF.matVLST, SURF.matCENTER, VISC.matFVLST, INPU.matROTORHUBGLOB, INPU.matROTORAXIS, SURF.matNTVLST] = fcnROTVEHICLE( SURF.matDVE, SURF.matVLST, SURF.matCENTER, INPU.valVEHICLES, SURF.vecDVEVEHICLE, INPU.matVEHORIG, VEHI.matVEHROT, VISC.matFVLST, VISC.vecFUSEVEHICLE, INPU.matROTORHUB, INPU.matROTORAXIS, VEHI.vecROTORVEH, SURF.matNTVLST, VISC.vecFDVEVEHICLE);

[ SURF.matUINF ] = fcnINITUINF( SURF.matCENTER, VEHI.matVEHUVW, VEHI.matVEHROT, SURF.vecDVEVEHICLE, ...
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

OUTP.DISTRIBUTIONS.vecDVEROTOR = SURF.vecDVEROTOR;
OUTP.DISTRIBUTIONS.matCENTER = SURF.matCENTER;
OUTP.DISTRIBUTIONS.vecDVEROTORBLADE = SURF.vecDVEROTORBLADE;

end
