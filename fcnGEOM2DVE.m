function [matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEROLL,...
    vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, matVLST0, matNTVLST0, matDVE, valNELE,...
    matADJE, vecDVESYM, vecDVETIP, vecDVESURFACE, vecDVELE, vecDVETE, vecDVEPANEL, matPANELTE,...
    valWINGS,vecDVEVEHICLE, vecDVEWING, vecDVEROTOR, vecDVEROTORBLADE, matSURFACETYPE, vecROTORVEH, ...
    matFUSEGEOM, matVEHUVW, matVEHROT, matVEHROTRATE, matCIRORIG, matUINF, vecDVETRI, vecN, vecM, valWSIZE, valWSIZETRI, vecPANELROTOR, vecQARM, cellAIRFOIL] = fcnGEOM2DVE(matGEOM, ...
    matVEHORIG, vecWINGTRI, vecWAKETRI, vecN, vecM, vecPANELWING,...
    vecSYM, vecSURFACEVEHICLE, vecROTOR, vecROTORBLADES, matROTORHUB, matROTORAXIS, matSECTIONFUSELAGE,...
    vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE, vecVEHVINF, vecVEHALPHA, vecVEHBETA, ...
    vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS, valVEHICLES, vecROTORRPM, vecPANELROTOR, cellAIRFOIL)

% tranlsate matGEOM to vehicle origin
matGEOM(:,1:3,:) = matGEOM(:,1:3,:)+permute(reshape(matVEHORIG(matGEOM(:,6,:),:)',3,2,[]),[2,1,3]);


vecPANELSURFACE = sort([vecPANELWING vecROTOR + max(vecPANELWING).*uint16((vecROTOR > 0))],2,'descend');
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
matDVE = uint16([]);
valNELE = 0;
matADJE = [];
vecDVESYM = uint8([]);
vecDVETIP = uint8([]);
vecDVESURFACE = uint8([]);
vecDVELE = uint8([]);
vecDVETE = uint8([]);
vecDVEPANEL = uint16([]);
matPANELTE = uint16([]);
vecDVETRI = [];

vecROTOR = uint8(vecROTOR);

valWSIZETRI = 0;
valWSIZE = 0;
tempM = [];
tempN = [];

for i = unique(vecPANELSURFACE,'stable')'
    
    panels = length(nonzeros(vecPANELSURFACE == i));
   
%     if isnan(vecWINGTRI(i))
        [tmatCENTER0, tvecDVEHVSPN, tvecDVEHVCRD, tvecDVELESWP, tvecDVEMCSWP, tvecDVETESWP, ...
            tvecDVEROLL, tvecDVEPITCH, tvecDVEYAW, tvecDVEAREA, tmatDVENORM, ...
            tmatVLST0, tmatNTVLST0, tmatDVE, tvalNELE, tmatADJE, ...
            tvecDVESYM, tvecDVETIP, tvecDVESURFACE, tvecDVELE, tvecDVETE, tvecDVEPANEL, tmatPANELTE] = ...
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
    
    vecDVESURFACE = [vecDVESURFACE; uint8(uint8(tvecDVESURFACE) + surfaceoffset)];
    vecDVEPANEL = [vecDVEPANEL; uint16(uint16(tvecDVEPANEL) + paneloffset)];
    
    vlstoffset = size(matVLST0,1);
    dveoffset = size(matDVE,1);
    matDVE = [matDVE; tmatDVE + vlstoffset];
    matVLST0 = [matVLST0; tmatVLST0];
    matNTVLST0 = [matNTVLST0; tmatNTVLST0];
  
    tmatADJE = [tmatADJE(:,1) + dveoffset tmatADJE(:,2) tmatADJE(:,3) + dveoffset tmatADJE(:,4)];
    matADJE = [matADJE; tmatADJE];
    temp = vecM(vecPANELSURFACE == i).*vecN(vecPANELSURFACE == i);
% 	tempM = [tempM; repmat(vecM(vecPANELSURFACE == i),[temp(1),1])];
%     tempN = [tempN; repmat(vecN(vecPANELSURFACE == i),[temp(1),1])];
end

% vecM = tempM;
% vecN = tempN;

% % Identifying which DVEs belong to which vehicle, as well as which type of lifting surface they belong to (wing or rotor)
vecDVEVEHICLE = vecSURFACEVEHICLE(vecDVESURFACE);
vecDVEWING = vecDVESURFACE;

vecDVEROTOR = vecROTOR(vecDVEPANEL); % Alton-Y
vecDVEROTORBLADE = vecDVEROTOR; % Current rotor DVEs are for Blade 1 (they are duplicated to Blade 2, 3, etc etc below)
idx_rotor = vecDVEROTOR>0; % Alton-Y
vecDVEWING(idx_rotor) = 0;

matSURFACETYPE = uint8(zeros(size(unique(vecDVESURFACE),1),2));
matSURFACETYPE(nonzeros(unique(vecDVEWING)),1) = nonzeros(unique(vecDVEWING));
matSURFACETYPE(nonzeros(unique(vecDVESURFACE(idx_rotor))),2) = nonzeros(unique(vecDVEROTOR));


% Identifying which ROTOR belongs to which vehicle.
vecROTORVEH = vecSURFACEVEHICLE(matSURFACETYPE(:,2)~=0);

% Duplicate Blades in a Rotor
[vecPANELROTOR, vecN, vecM, matVLST0, matCENTER0, matDVE, matADJE, vecDVEVEHICLE, ...
    vecDVEWING, vecDVEROTOR, matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, ...
    vecDVETIP, vecDVELE, vecDVETE, vecDVEROTORBLADE, vecDVESYM, ...
    valNELE, matNTVLST0, cellAIRFOIL] = fcnDUPBLADE( vecROTORVEH, vecDVEROTOR, ...
    matVLST0, matCENTER0, matDVE, matADJE, vecROTORBLADES, ...
    valNELE, matROTORHUB, matVEHORIG, vecDVEVEHICLE, vecDVEWING, ...
    matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, vecDVETIP, vecDVELE, ...
    vecDVETE, vecDVEROTORBLADE, vecDVESYM, matROTORAXIS, matNTVLST0, ...
    vecM, vecN, vecPANELROTOR, cellAIRFOIL);

matFUSEGEOM = fcnCREATEFUSE(matSECTIONFUSELAGE, vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE);


% flap = 40;
% idx = find(vecDVETE > 0 & vecDVEWING > 0);
% for i = idx'
% u = matVLST0(matDVE(i,2),:) - matVLST0(matDVE(i,1),:);
% uo = (matVLST0(matDVE(i,2),:) + matVLST0(matDVE(i,1),:))./2;
% 
% R = columbia_rotation(u, -flap);
% matVLST0(matDVE(i,3),:) = (matVLST0(matDVE(i,3),:) - uo)*R + uo;
% matVLST0(matDVE(i,4),:) = (matVLST0(matDVE(i,4),:) - uo)*R + uo;
% end



[ matVEHUVW, matVEHROT, matVEHROTRATE, matCIRORIG, vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, matVEHORIG, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS );
[matVLST0, matCENTER0, matFUSEGEOM, matROTORHUBGLOB, matROTORAXIS, matNTVLST0] = fcnROTVEHICLE( matDVE, matVLST0, matCENTER0, valVEHICLES, vecDVEVEHICLE, matVEHORIG, matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS, matROTORHUB, matROTORAXIS, vecROTORVEH, matNTVLST0);

[ matUINF ] = fcnINITUINF( matCENTER0, matVEHUVW, matVEHROT, vecDVEVEHICLE, ...
    vecDVEROTOR, vecROTORVEH, matVEHORIG, matROTORHUBGLOB, matROTORAXIS, vecROTORRPM );


% update DVE params after vehicle rotation
[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, matCENTER0 ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST0);

valWSIZE = length(nonzeros(vecDVETE));

% Compute torque arm length for rotor power calculations
vecQARM = zeros(valNELE,3);
if max(vecDVEROTOR)>0
    vecQARM(vecDVEROTOR>0,:) = matCENTER0(vecDVEROTOR>0,:) - matROTORHUB(vecDVEROTOR(vecDVEROTOR>0),:);    
    vecQARM = sqrt(sum(vecQARM.^2,2));
end

end

function [R] = columbia_rotation(u,theta)
%ROTATION Summary of this function goes here
%   Detailed explanation goes here
ux = u(1);
uy = u(2);
uz = u(3);

cost = cosd(theta);
sint = sind(theta);

R = [cost + ux.^2*(1 - cost) ux*uy*(1 - cost) - uz*sint ux*uz*(1 - cost) + uy*sint; ...
    uy*ux*(1 - cost) + uz*sint cost + uy.^2*(1 - cost) uy*uz*(1 - cost) - ux*sint; ...
    uz*ux*(1 - cost) - uy*sint uz*uy*(1 - cost) + ux*sint cost + uz.^2*(1 - cost); ...
    ];

end
