function [FLAG, COND, VISC, INPU, VEHI, SURF] = fcnXMLREAD(filename, VAP_IN)

% clc
% clear
%
% filename = 'inputs/XMLtest.vap';

% OUTPUT

% vecROTORS - Rows are wing number, tells which vehicle has how many rotors
% COND.vecROTORRPM - Rows are rotor number, tells the rotor RPM
% INPU.vecROTDIAM - Rows are rotor number, tells the rotor diameter
% vecROTORHUB - Rows are rotor number, column is xyz in vehicle local. Tells rotor hub location
% INPU.matROTORAXIS - Rows are rotor number, column is xyz in vehicle local. Unit vector of rotor plane normal
% vecROTORM - Rows are rotor number, tells how many chordwise elements on the rotor blade
% vecROTOR - Rows are panel number, tells which rotor the panel belongs to

% MISC.matSURFACETYPE - surface x 2 matrix defining which surface is a wing or rotor, where 1 in column 1 is wing, 1 in column 2 is rotor

% VEHI.vecSURFACEVEHICLE - Rows are wing number, tells which vehicle each wing belongs to
% vecWING - Rows are panel number, tells which wing each panel belongs to
% INPU.valPANELS - Total number of panels
% vecPANELS - Rows are wing number, tells how many panels are on each wing
% vecWINGM - Rows are wing number, tells us the number of chordwise lifting lines on each wing
% vecWINGCMAC - Rows are wing number, tells us the mean aerodynamic chord of each wing
% vecWINGSPAN - Rows are wing number, tell us the span of each wing (tip to tip)
% INPU.vecSYM - Rows are panel number, tells us whether edge 1 or 2 of the panel has symmetry boundary condition
% INPU.vecN - Rows are panel number, tells us the number of spanwise DVEs for each panel
% INPU.vecM - Rows are panel number, tells us the number of chordwise DVEs for each panel (same as vecWINGM but on the DVE level)
% vecWINGINCID - Rows are wing number, tells us the incidence angle of the wing
% matSECTIONS - number of sections x 5 matrix of [x y z chord twist] for each section

%%
inp = fcnXML2STRUCT(filename);
VAP = inp.VAP;

%% Settings
if strcmpi(VAP.settings.relax.Text, 'true') FLAG.RELAX = 1; else FLAG.RELAX = 0; end
if strcmpi(VAP.settings.steady.Text, 'true') FLAG.STEADY = 1; else FLAG.STEADY = 0; end
try
    if strcmpi(VAP.settings.stiff_wing.Text, 'true') FLAG.STIFFWING = 1;
    else FLAG.STIFFWING = 2;
    end
catch
    FLAG.STIFFWING = 1;
end

try if strcmpi(VAP.settings.fixed_lift.Text, 'true') FLAG.FIXEDLIFT = 1; end; catch FLAG.FIXEDLIFT = 0; end
try FLAG.GUSTMODE = str2double(VAP.settings.gust_mode.Text); catch; FLAG.GUSTMODE = 0; end
if strcmpi(VAP.settings.gliding.Text, 'true') FLAG.GLIDING = 1; else; FLAG.GLIDING = 0; end

COND.valMAXTIME = floor(str2double(VAP.settings.maxtime.Text));
COND.valDELTIME = str2double(VAP.settings.delta_time.Text);
try COND.valSTARTFORCES = floor(str2double(VAP.settings.start_forces.Text)); catch; COND.valSTARTFORCES = 0; end

%% Conditions
COND.valDENSITY = str2double(VAP.conditions.density.Text);
VISC.valKINV = str2double(VAP.conditions.kin_viscosity.Text);
try COND.valGUSTAMP = str2double(VAP.conditions.gust_amplitude.Text); catch; COND.valGUSTAMP = 0; end
try COND.valGUSTL = str2double(VAP.conditions.gust_length.Text); catch; COND.valGUSTL = 0; end
try COND.valGUSTSTART = str2double(VAP.conditions.gust_start.Text); catch; COND.valGUSTSTART = 0; end

%% Vehicles
INPU.valVEHICLES = max(size(VAP.vehicle));

INPU.matVEHORIG = nan(INPU.valVEHICLES,3);
INPU.vecVEHCG = nan(INPU.valVEHICLES,3);
COND.vecVEHVINF = nan(INPU.valVEHICLES,1);
COND.vecVEHALPHA = nan(INPU.valVEHICLES,1);
COND.vecVEHPITCH = nan(INPU.valVEHICLES,1);
COND.vecVEHBETA = nan(INPU.valVEHICLES,1);
COND.vecVEHROLL = nan(INPU.valVEHICLES,1);
COND.vecVEHFPA = nan(INPU.valVEHICLES,1);
COND.vecVEHTRK = nan(INPU.valVEHICLES,1);
COND.vecVEHWEIGHT = nan(INPU.valVEHICLES,1);

VEHI.vecVEHRADIUS = nan(INPU.valVEHICLES,1);

vecWINGS = nan(INPU.valVEHICLES,1);
vecSTRUCTURE = nan(INPU.valVEHICLES,1);
vecROTORS = nan(INPU.valVEHICLES,1);
vecFUSELAGES = nan(INPU.valVEHICLES,1);

COND.vecROTORRPM = [];
INPU.vecROTDIAM = [];
INPU.matROTORHUB = [];
INPU.matROTORAXIS = [];
INPU.vecROTORBLADES = [];
vecROTORM = [];
COND.vecCOLLECTIVE = [];

VISC.vecINTERF = 0;

vecWINGINCID = [];
FLAG.vecTRIMABLE = [];
vecWINGM = [];
SURF.matWINGORIG = [];
vecPANELS = [];
cellAIRFOILtemp = {};
VISC.cellAIRFOIL = {};
vecSECTIONS = [];
matSECTIONS = [];
vecSECTIONPANEL = [];

COND.vecSURFTRI = [];

k = 1;
kk = 1;
kkk = 1;
kkkk = 1;

o = 1;
p = 1;

q = 1;

for i = 1:INPU.valVEHICLES
    
    try veh = VAP.vehicle{1,i}; catch; veh = VAP.vehicle; end
    
    INPU.matVEHORIG(i,:) = [str2double(veh.global_x.Text) str2double(veh.global_y.Text) str2double(veh.global_z.Text)];
    try COND.vecVEHVINF(i,1) = str2double(veh.speed.Text); catch; COND.vecVEHVINF(i,1) = nan; end
    try COND.vecVEHWEIGHT(i,1) = str2double(veh.weight.Text); catch; COND.vecVEHWEIGHT(i,1) = nan; end
    COND.vecVEHALPHA(i,1) = str2double(veh.alpha.Text);
    COND.vecVEHBETA(i,1) = str2double(veh.beta.Text);
    COND.vecVEHROLL(i,1) = str2double(veh.roll.Text);
    COND.vecVEHFPA(i,1) = str2double(veh.fpa.Text);
    COND.vecVEHTRK(i,1) = str2double(veh.track.Text);
    try INPU.vecVEHCG(i,:) = [str2double(veh.vehicle_CG.x.Text) str2double(veh.vehicle_CG.y.Text) str2double(veh.vehicle_CG.z.Text)]; end
    
    try VEHI.vecVEHRADIUS(i,1) = str2double(veh.radius.Text); catch VEHI.vecVEHRADIUS(i,1) = nan; end
    
    try vecWINGS(i,1) = max(size(veh.wing)); catch; vecWINGS(i,1) = 0; end
    try vecSTRUCTURE(i,1) = max(size(veh.structure)); catch; vecSTRUCTURE(i,1) = 0; end
    try vecROTORS(i,1) = max(size(veh.rotor)); catch; vecROTORS(i,1) = 0; end
    try vecFUSELAGES(i,1) = max(size(veh.fuselage)); catch; vecFUSELAGES(i,1) = 0; end
    try vecPAYLOADS(i,1) = max(size(veh.payload)); catch; vecPAYLOADS(i,1) = 0; end
    try vecPROPS(i,1) = max(size(veh.propeller)); catch; vecPROPS(i,1) = 0; end
    
    INPU.vecAREA(i) = str2double(veh.ref_area.Text);
    INPU.vecSPAN(i) = str2double(veh.ref_span.Text);
    INPU.vecCMAC(i) = str2double(veh.ref_cmac.Text);
    
    try VISC.vecINTERF = str2double(veh.interference_drag.Text); end
    
    %% Loading Wings
    for j = 1:vecWINGS(i)
        
        try win = veh.wing{1,j}; catch; win = veh.wing; end
        
        COND.vecSURFTRI(k,1) = false;
        try if strcmpi(win.triangular_elements.Text, 'true') COND.vecSURFTRI(k) = true; end; end;
        try if strcmpi(win.symmetry.Text, 'true') sym = true; else sym = false; end; catch sym = false; end;
        
        try if strcmpi(win.type.Text, 'main wing') INPU.vecWINGTYPE(k) = 1; elseif strcmpi(win.type.Text, 'hstab') INPU.vecWINGTYPE(k) = 2; elseif strcmpi(win.type.Text, 'vstab') INPU.vecWINGTYPE(k) = 3; end; catch INPU.vecWINGTYPE(k) = 1; end;
        
        vecWINGINCID(k) = str2double(win.incidence.Text);
        if strcmpi(win.trimable.Text, 'true') FLAG.vecTRIMABLE(j,i) = 1; else FLAG.vecTRIMABLE(j,i) = 0; end
        try if strcmpi(win.flexible.Text, 'true') FLAG.vecFLEXIBLE(k,i) = 1; else FLAG.vecFLEXIBLE(k,i) = 0; end; catch; FLAG.vecFLEXIBLE(k,i) = 0; end
        
        vecWINGM(k,1) = str2double(win.chordwise_elements.Text);
        
        try matWINGORIG(k,:) = [str2double(win.vehicle_x.Text) str2double(win.vehicle_y.Text) str2double(win.vehicle_z.Text)];
        catch; matWINGORIG(k,:) = [0 0 0]; end
        
        try SURF.matTRIMORIG(k,:) = [str2double(win.trimorigin.x.Text) str2double(win.trimorigin.y.Text) str2double(win.trimorigin.z.Text)];
        catch; SURF.matTRIMORIG(k,:) = [0 0 0]; end
        
        try VEHI.vecWINGMASS(k,:) = str2double(win.mass.Text);
        catch; VEHI.vecWINGMASS(k,:) = 0; end
        
        try VEHI.vecWINGCG(k,:) = [str2double(win.CG.x.Text) str2double(win.CG.y.Text) str2double(win.CG.z.Text)];
        catch; VEHI.vecWINGCG(k,:) = [0 0 0]; end
        
        vecPANELS(k,1) = max(size(win.panel));
        
        for m = 1:vecPANELS(k,1)
            
            try pan = win.panel{1,m}; catch; pan = win.panel; end
            
            vecSYMtemp(kk,1) = sym;
            
            try cellAIRFOILtemp{kk} = pan.strip_airfoil.Text; catch cellAIRFOILtemp{kk} = nan; end
            vecNtemp(kk,1) = floor(str2double(pan.spanwise_elements.Text));
            vecMtemp(kk,1) = str2double(win.chordwise_elements.Text);
            try cellNSPACEtemp{kk,1} = pan.spanwise_spacing.Text; catch cellNSPACEtemp{kk,1} = 'NORMAL'; end
            try cellMSPACEtemp{kk,1} = win.chordwise_spacing.Text; catch cellMSPACEtemp{kk,1} = 'NORMAL'; end   
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                                
                if vecSECTIONS(kk,1) > 1
                    sec = pan.section{1,n};
                else
                    sec = pan.section(1,n);
                end
                
                matSECTIONS(kkk,:) = [str2double(sec.wing_x.Text) + matWINGORIG(k,1) str2double(sec.wing_y.Text) + matWINGORIG(k,2) str2double(sec.wing_z.Text) + matWINGORIG(k,3) str2double(sec.chord.Text) vecWINGINCID(k)+str2double(sec.twist.Text) i];
                try cellCAIRFOIL{kkk,:} = sec.camber_airfoil.Text; catch cellCAIRFOIL{kkk,:} = nan; end
                
                vecSECTIONPANEL(kkk,1) = kk;
                
                kkk = kkk + 1;
            end
            
            INPU.vecPANELWING(kk,1) = o;
            panel_rotors(kk,1) = 0;
            
            kk = kk + 1;
        end
        
        o = o + 1;
        
        VEHI.vecSURFACEVEHICLE(k,1) = i;
        k = k + 1;
    end
    
    %% Loading Wing Structure
    
    for j = 1:vecSTRUCTURE(i)
        
        try struct = VAP.vehicle.structure{1,i}; catch; struct = VAP.vehicle.structure; end
        
        COND.valSDELTIME = str2double(struct.conditions.valSDELTIME.Text);
        INPU.valNSELE = str2double(struct.conditions.valNSELE.Text);
        COND.valSTIFFSTEPS = str2double(struct.conditions.valSTIFFSTEPS.Text);
        COND.valSTAGGERSTEPS = str2double(struct.conditions.valSTAGGERSTEPS.Text);
        
        INPU.vecEIxCOEFF(1) = str2double(struct.properties.stiffness.A_vecEIx.Text);
        INPU.vecEIxCOEFF(2) = str2double(struct.properties.stiffness.B_vecEIx.Text);
        INPU.vecEIxCOEFF(3) = str2double(struct.properties.stiffness.C_vecEIx.Text);
        
        INPU.vecGJtCOEFF(1) = str2double(struct.properties.stiffness.A_vecGJt.Text);
        INPU.vecGJtCOEFF(2) = str2double(struct.properties.stiffness.B_vecGJt.Text);
        INPU.vecGJtCOEFF(3) = str2double(struct.properties.stiffness.C_vecGJt.Text);
        
        INPU.vecEACOEFF(1) = str2double(struct.properties.geometry.A_vecEA.Text);
        INPU.vecEACOEFF(2) = str2double(struct.properties.geometry.B_vecEA.Text);
        INPU.vecEACOEFF(3) = str2double(struct.properties.geometry.C_vecEA.Text);
        
        INPU.vecCGCOEFF(1) = str2double(struct.properties.geometry.A_vecCG.Text);
        INPU.vecCGCOEFF(2) = str2double(struct.properties.geometry.B_vecCG.Text);
        INPU.vecCGCOEFF(3) = str2double(struct.properties.geometry.C_vecCG.Text);
        
        INPU.vecJTCOEFF(1) = str2double(struct.properties.mass.A_vecJt.Text);
        INPU.vecJTCOEFF(2) = str2double(struct.properties.mass.B_vecJt.Text);
        INPU.vecJTCOEFF(3) = str2double(struct.properties.mass.C_vecJt.Text);
        
        INPU.vecLMCOEFF(1) = str2double(struct.properties.mass.A_vecLM.Text);
        INPU.vecLMCOEFF(2) = str2double(struct.properties.mass.B_vecLM.Text);
        INPU.vecLMCOEFF(3) = str2double(struct.properties.mass.C_vecLM.Text);
        
    end
    
    %% Loading Vehicle Fuselage
    
    for j = 1:vecFUSELAGES(i)
        
        try fuse = VAP.vehicle.fuselage{1,i}; catch; fuse = VAP.vehicle.fuselage; end
        
        try VEHI.vecFUSEMASS(j,:) = str2double(fuse.mass.Text);
        catch; VEHI.vecFUSEMASS(j,:) = 0; end
        
        try VEHI.vecFUSEL(j,:) = str2double(fuse.length.Text);
        catch; VEHI.vecFUSEMASS(j,:) = 0; end
        try VEHI.valNFELE(j,:) = str2double(fuse.valNFELE.Text);
        catch; VEHI.valNFELE(j,:) = 0; end
        
        try VEHI.vecFUSELOC(j,:) = [str2double(fuse.start_loc.x.Text) str2double(fuse.start_loc.y.Text) str2double(fuse.start_loc.z.Text)];
        catch; VEHI.vecFUSELOC(j,:) = [0 0 0]; end
        
    end
    
    %% Loading Vehicle Payload
    
    for j = 1:vecPAYLOADS(i)
        
        try payl = VAP.vehicle.payload{1,i}; catch; payl = VAP.vehicle.payload; end
        
        try VEHI.vecPAYLMASS(j,:) = str2double(payl.mass.Text);
        catch; VEHI.vecPAYLMASS(j,:) = 0; end
        
        try VEHI.vecPAYLCG(j,:) = [str2double(payl.CG.x.Text) str2double(payl.CG.y.Text) str2double(payl.CG.z.Text)];
        catch; VEHI.vecPAYLCG(j,:) = [0 0 0]; end
        
    end
    
    %% Loading Propeller (if not being modelled in the aerodynamic model)
    
    for j = 1:vecPROPS(i)
        
        try prop = VAP.vehicle.propeller{1,i}; catch; prop = VAP.vehicle.propeller; end
        
        try VEHI.vecPROPLOC(j,:) = [str2double(prop.location.x.Text) str2double(prop.location.y.Text) str2double(prop.location.z.Text)];
        catch; VEHI.vecPROPLOC(j,:) = []; end
        
        try VEHI.vecPROPDIR(j,:) = [str2double(prop.direction.x.Text) str2double(prop.direction.y.Text) str2double(prop.direction.z.Text)];
        catch; VEHI.vecPROPDIR(j,:) = []; end
        
    end
    
    %% Loading Rotors
    for j = 1:vecROTORS(i)
        
        try rot = veh.rotor{1,j}; catch; rot = veh.rotor; end
        COND.vecSURFTRI(k,1) = false;
        try if strcmpi(rot.triangular_elements.Text, 'true') COND.vecSURFTRI(k) = true; end; end;
        try if strcmpi(rot.symmetry.Text, 'true') sym = true; else sym = false; end; catch sym = false; end;
        
        try if strcmpi(win.type.Text, 'rotor') INPU.vecWINGTYPE(k) = 4; end; catch INPU.vecWINGTYPE(k) = 4; end;
        
        COND.vecROTORRPM(p,1) = str2double(rot.rpm.Text);
        INPU.vecROTDIAM(p,1) = str2double(rot.ref_diam.Text);
        
        INPU.matROTORHUB(p,:) = [str2double(rot.veh_x_hub.Text) str2double(rot.veh_y_hub.Text) str2double(rot.veh_z_hub.Text)];
        INPU.matROTORAXIS(p,:) = [str2double(rot.veh_x_axis.Text) str2double(rot.veh_y_axis.Text) str2double(rot.veh_z_axis.Text)];
        
        INPU.vecROTORBLADES(p,:) = floor(str2double(rot.blades.Text));
        
        vecROTORM(p,1) = floor(str2double(rot.chordwise_elements.Text));
        
        try COND.vecCOLLECTIVE(p,1) = str2double(rot.collective.Text); catch; COND.vecCOLLECTIVE(p,1) = 0; end;
        
        vecPANELS(k,1) = max(size(rot.panel));
        
        flip = 1;
        try if strcmpi(rot.rotation_direction.Text,'CW'); flip = -1; end; end
        COND.vecROTORRPM(p,1) = COND.vecROTORRPM(p,1)*flip;
        for m = 1:vecPANELS(k,1)
            
            try pan = rot.panel{1,m}; catch; pan = rot.panel; end
            
            vecSYMtemp(kk,1) = sym;
            
            try cellAIRFOILtemp{kk} = pan.strip_airfoil.Text; catch cellAIRFOILtemp{kk} = nan; end
            vecNtemp(kk,1) = floor(str2double(pan.spanwise_elements.Text));
            vecMtemp(kk,1) = str2double(rot.chordwise_elements.Text);
            try cellNSPACEtemp{kk,1} = pan.spanwise_spacing.Text; catch cellNSPACEtemp{kk,1} = 'NORMAL'; end
            try cellMSPACEtemp{kk,1} = rot.chordwise_spacing.Text; catch cellMSPACEtemp{kk,1} = 'NORMAL'; end   
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                sec = pan.section{1,n};
                
                matSECTIONS(kkk,:) = [str2double(sec.rotor_x.Text) + INPU.matROTORHUB(p,1) flip*str2double(sec.rotor_y.Text) + INPU.matROTORHUB(p,2) str2double(sec.rotor_z.Text) + INPU.matROTORHUB(p,3) str2double(sec.chord.Text) str2double(sec.twist.Text) i];
                try cellCAIRFOIL{kkk,:} = sec.camber_airfoil.Text; catch cellCAIRFOIL{kkk,:} = nan; end
                vecSECTIONPANEL(kkk,1) = kk;
                
                kkk = kkk + 1;
            end
            
            INPU.vecPANELWING(kk,1) = 0;
            panel_rotors(kk,1) = p;
            
            kk = kk + 1;
        end
        
        p = p + 1;
        
        VEHI.vecSURFACEVEHICLE(k,1) = i;
        k = k + 1;
    end
end

INPU.valPANELS = sum(vecPANELS);

%% Reorganizing data for output
k = 1;
for i = 1:INPU.valPANELS
    sections = matSECTIONS(vecSECTIONPANEL == i,:);
    len = size(sections,1);

    for j = 1:len - 1
        INPU.matGEOM(1,:,k) = sections(j,:);
        INPU.matGEOM(2,:,k) = sections(j+1,:);
        INPU.cellCAIRFOIL{k,1} = cellCAIRFOIL{j};
        INPU.cellCAIRFOIL{k,2} = cellCAIRFOIL{j+1};
        INPU.vecPANELWING(k,1) = INPU.vecPANELWING(i);
        INPU.vecPANELROTOR(k,1) = panel_rotors(i);

        INPU.vecSYM(k,1) = vecSYMtemp(i);           
        INPU.vecN(k,1) = vecNtemp(i);
        INPU.vecM(k,1) = vecMtemp(i);
        INPU.cellMSPACE{k,1} = cellMSPACEtemp{i,1};
        INPU.cellNSPACE{k,1} = cellNSPACEtemp{i,1};

        if ~isempty(cellAIRFOILtemp)
            VISC.cellAIRFOIL{k,1} = cellAIRFOILtemp{i};
        end

        k = k + 1;
    end     
end

INPU.valPANELS = size(INPU.matGEOM,3);

%% Overwriting default values with user-specified

if ~isempty(VAP_IN)
    in_names = fieldnames(VAP_IN);
    for i = 1:length(in_names)
        if isfield(INPU, in_names{i})
            INPU.(in_names{i}) = VAP_IN.(in_names{i});
        elseif isfield(COND, in_names{i})
            COND.(in_names{i}) = VAP_IN.(in_names{i});
        elseif isfield(VEHI, in_names{i})
            VEHI.(in_names{i}) = VAP_IN.(in_names{i});
        elseif isfield(VISC, in_names{i})
            VISC.(in_names{i}) = VAP_IN.(in_names{i});
        elseif isfield(FLAG, in_names{i})
            FLAG.(in_names{i}) = VAP_IN.(in_names{i});
        end
    end
end

if any(isnan(COND.vecVEHVINF))
%     disp('One or more vehicle velocities was not read in - setting to unity and enabling fixed-lift analysis');
    COND.vecVEHVINF(isnan(COND.vecVEHVINF)) = ones(sum(isnan(COND.vecVEHVINF)),1);
    FLAG.FIXEDLIFT = 1;
else
    FLAG.FIXEDLIFT = 0;
end

end








