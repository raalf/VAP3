function [matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEROLL,...
    vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, matVLST0, matNTVLST0, matNPVLST0, matDVE, valNELE,...
    matADJE, vecDVESYM, vecDVETIP, vecDVESURFACE, vecDVELE, vecDVETE, vecDVEPANEL, matPANELTE,...
    valWINGS,vecDVEVEHICLE, vecDVEWING, vecDVEROTOR, vecDVEROTORBLADE, matSURFACETYPE, vecROTORVEH, ...
    matFUSEGEOM, matVEHUVW, matVEHROT, matVEHROTRATE, matCIRORIG, vecVEHPITCH, vecVEHYAW,...
    matROTORHUBGLOB, matUINF, vecDVETRI, vecN, vecM, valWSIZE, valWSIZETRI] = fcnGEOM2DVE(matGEOM, ...
    matVEHORIG, vecWINGTRI, vecWAKETRI, vecN, vecM, vecPANELWING,...
    vecSYM, vecSURFACEVEHICLE, vecROTOR, vecROTORBLADES, matROTORHUB, matROTORAXIS, matSECTIONFUSELAGE,...
    vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE, vecVEHVINF, vecVEHALPHA, vecVEHBETA, ...
    vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS, valVEHICLES, vecROTORRPM)

% tranlsate matGEOM to vehicle origin
matGEOM(:,1:3,:) = matGEOM(:,1:3,:)+permute(reshape(matVEHORIG(matGEOM(:,6,:),:)',3,2,[]),[2,1,3]);


vecPANELSURFACE = sort([vecPANELWING vecROTOR],2,'descend');
vecPANELSURFACE = vecPANELSURFACE(:,1);

valWINGS = max(vecPANELSURFACE);

matCENTER0 = [];
vecDVEHVSPN = [];
vecDVEHVCRD = [];
vecDVELESWP = [];
vecDVEMCSWP = [];
vecDVETESWP = [];
vecDVEROLL = [];
vecDVEPITCH = [];
vecDVEYAW = [];
vecDVEAREA = [];
matDVENORM = [];
matVLST0 = [];
matNTVLST0 = [];
matNPVLST0 = [];
matDVE = [];
valNELE = 0;
matADJE = [];
vecDVESYM = [];
vecDVETIP = [];
vecDVESURFACE = [];
vecDVELE = [];
vecDVETE = [];
vecDVEPANEL = [];
matPANELTE = [];
vecDVETRI = [];

valWSIZETRI = 0;
valWSIZE = 0;

for i = 1:valWINGS
    
    panels = length(nonzeros(vecPANELSURFACE == i));
    
%     if isnan(vecWINGTRI(i))
        [tmatCENTER0, tvecDVEHVSPN, tvecDVEHVCRD, tvecDVELESWP, tvecDVEMCSWP, tvecDVETESWP, ...
            tvecDVEROLL, tvecDVEPITCH, tvecDVEYAW, tvecDVEAREA, tmatDVENORM, ...
            tmatVLST0, tmatNTVLST0, tmatDVE, tvalNELE, tmatADJE, ...
            tvecDVESYM, tvecDVETIP, tvecDVESURFACE, tvecDVELE, tvecDVETE, tvecDVEPANEL, tmatNPVLST0, tmatPANELTE] = ...
            fcnGENERATEDVES(panels, matGEOM(:,:,(vecPANELSURFACE == i)), vecSYM(vecPANELSURFACE == i), vecN(vecPANELSURFACE == i), vecM(vecPANELSURFACE == i));
        
        vecDVETRI = [vecDVETRI; zeros(tvalNELE,1)];
        
%     else
%         
%         [tmatCENTER0, tvecDVEHVSPN, tvecDVEHVCRD, tvecDVELESWP, tvecDVEMCSWP, tvecDVETESWP, ...
%             tvecDVEROLL, tvecDVEPITCH, tvecDVEYAW, tvecDVEAREA, tmatDVENORM, ...
%             tmatVLST0, tmatNTVLST0, tmatDVE, tvalNELE, tmatADJE, ...
%             tvecDVESYM, tvecDVETIP, tvecDVESURFACE, tvecDVELE, tvecDVETE, ...
%             tvecDVEPANEL, tmatNPVLST0, vecM(vecPANELSURFACE == i,1), vecN(vecPANELSURFACE == i), tmatPANELTE] = ...
%             fcnGENERATEDVESTRI(panels, matGEOM(:,:,(vecPANELSURFACE == i)), vecSYM(vecPANELSURFACE == i), vecN(vecPANELSURFACE == i), vecM(vecPANELSURFACE == i));
%         
%         vecDVETRI = [vecDVETRI; ones(tvalNELE,1)];
        
        
%     end
    
%     if isnan(vecWAKETRI(i))

%     else
%         valWSIZETRI = valWSIZETRI + sum(nonzeros(tvecDVETE > 0))*2;
%     end
    
    valNELE = valNELE + tvalNELE;
    matCENTER0 = [matCENTER0; tmatCENTER0];
    vecDVEHVSPN = [vecDVEHVSPN; tvecDVEHVSPN];
    vecDVEHVCRD = [vecDVEHVCRD; tvecDVEHVCRD];
    vecDVELESWP = [vecDVELESWP; tvecDVELESWP];
    vecDVEMCSWP = [vecDVEMCSWP; tvecDVEMCSWP];
    vecDVETESWP = [vecDVETESWP; tvecDVETESWP];
    vecDVEROLL = [vecDVEROLL; tvecDVEROLL];
    vecDVEPITCH = [vecDVEPITCH; tvecDVEPITCH];
    vecDVEYAW = [vecDVEYAW; tvecDVEYAW];
    vecDVEAREA = [vecDVEAREA; tvecDVEAREA];
    matDVENORM = [matDVENORM; tmatDVENORM];
    vecDVESYM = [vecDVESYM; tvecDVESYM];
    vecDVETIP = [vecDVETIP; tvecDVETIP];
    vecDVELE = [vecDVELE; tvecDVELE];
    vecDVETE = [vecDVETE; tvecDVETE];
    matPANELTE = [matPANELTE; tmatPANELTE];
    
    if i == 1; surfaceoffset = 0; paneloffset = 0;
    else; surfaceoffset = max(vecDVESURFACE); paneloffset = max(vecDVEPANEL);
    end
    
    vecDVESURFACE = [vecDVESURFACE; tvecDVESURFACE + surfaceoffset];
    vecDVEPANEL = [vecDVEPANEL; tvecDVEPANEL + paneloffset];
    
    vlstoffset = size(matVLST0,1);
    dveoffset = size(matDVE,1);
    matDVE = [matDVE; tmatDVE + vlstoffset];
    matVLST0 = [matVLST0; tmatVLST0];
    matNTVLST0 = [matNTVLST0; tmatNTVLST0];
    matNPVLST0 = [matNPVLST0; tmatNPVLST0];
    
    tmatADJE = [tmatADJE(:,1) + dveoffset tmatADJE(:,2) tmatADJE(:,3) + dveoffset tmatADJE(:,4)];
    matADJE = [matADJE; tmatADJE];
    
end

% % Identifying which DVEs belong to which vehicle, as well as which type of lifting surface they belong to (wing or rotor)
vecDVEVEHICLE = vecSURFACEVEHICLE(vecDVESURFACE);
vecDVEWING = vecDVESURFACE;

vecDVEROTOR = vecROTOR(vecDVEPANEL); % Alton-Y
vecDVEROTORBLADE = vecDVEROTOR; % Current rotor DVEs are for Blade 1 (they are duplicated to Blade 2, 3, etc etc below)
idx_rotor = vecDVEROTOR>0; % Alton-Y
vecDVEWING(idx_rotor) = 0;

matSURFACETYPE = zeros(size(unique(vecDVESURFACE),1),2);
matSURFACETYPE(nonzeros(unique(vecDVEWING)),1) = nonzeros(unique(vecDVEWING));
matSURFACETYPE(nonzeros(unique(vecDVESURFACE(idx_rotor))),2) = nonzeros(unique(vecDVEROTOR));


% Identifying which ROTOR belongs to which vehicle.
vecROTORVEH = vecSURFACEVEHICLE(matSURFACETYPE(:,2)~=0);

% Duplicate Blades in a Rotor
[ matVLST0, matNPVLST0, matCENTER0, matDVE, matADJE, vecDVEVEHICLE, ...
    vecDVEWING, vecDVEROTOR, matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, ...
    vecDVETIP, vecDVELE, vecDVETE, vecDVEROTORBLADE, vecDVESYM, ...
    valNELE, matNTVLST0] = fcnDUPBLADE( vecROTORVEH, vecDVEROTOR, ...
    matVLST0, matCENTER0, matNPVLST0, matDVE, matADJE, vecROTORBLADES, ...
    valNELE, matROTORHUB, matVEHORIG, vecDVEVEHICLE, vecDVEWING, ...
    matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, vecDVETIP, vecDVELE, ...
    vecDVETE, vecDVEROTORBLADE, vecDVESYM, matROTORAXIS, matNTVLST0);


matFUSEGEOM = fcnCREATEFUSE(matSECTIONFUSELAGE, vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE);

% DOES THIS HAVE TO GO BEFORE DUPBLADE??? MATNPVLST IS GETTING FUCKED BY A TRANSFORMATION FROM LOCAL ROTOR TO LCOAL VEHICL WHICH WE ASSUME IS GLOBAL???
[ matVEHUVW, matVEHROT, matVEHROTRATE, matCIRORIG, vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, matVEHORIG, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS );
[matVLST0, matCENTER0, matFUSEGEOM, matROTORHUBGLOB, matROTORAXIS, matNPVLST0, matNTVLST0] = fcnROTVEHICLE( matDVE, matVLST0, matCENTER0, valVEHICLES, vecDVEVEHICLE, matVEHORIG, matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS, matROTORHUB, matROTORAXIS, vecROTORVEH, matNPVLST0, matNTVLST0);

[ matUINF ] = fcnINITUINF( matCENTER0, matVEHUVW, matVEHROT, vecDVEVEHICLE, ...
    vecDVEROTOR, vecROTORVEH, matVEHORIG, matROTORHUBGLOB, matROTORAXIS, vecROTORRPM );


% update DVE params after vehicle rotation
[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, ~ ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST0);

valWSIZE = length(nonzeros(vecDVETE));

end

