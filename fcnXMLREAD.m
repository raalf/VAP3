function [FLAG, COND, VISC, INPU, VEHI] = fcnXMLREAD(filename, VAP_IN)

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
if strcmpi(VAP.settings.flagRELAX.Text, 'true') FLAG.RELAX = 1; else FLAG.RELAX = 0; end
if strcmpi(VAP.settings.flagSTEADY.Text, 'true') FLAG.STEADY = 1; else FLAG.STEADY = 0; end
try
    if strcmpi(VAP.settings.flagSTIFFWING.Text, 'true') FLAG.STIFFWING = 1;
    else FLAG.STIFFWING = 2;
    end
catch
    FLAG.STIFFWING = 1;
end
% if strcmpi(VAP.settings.FLAG.TRI.Text, 'true') FLAG.TRI = 1; else FLAG.TRI = 0; end

try if strcmpi(VAP.settings.flagFIXEDLIFT.Text, 'true') FLAG.FIXEDLIFT = 1; end; catch FLAG.FIXEDLIFT = 0; end
try FLAG.GUSTMODE = str2double(VAP.conditions.flagGUSTMODE.Text); catch; FLAG.GUSTMODE = 0; end

COND.valMAXTIME = floor(str2double(VAP.settings.valMAXTIME.Text));
COND.valMINTIME = floor(str2double(VAP.settings.valMINTIME.Text));
COND.valDELTIME = str2double(VAP.settings.valDELTIME.Text);
COND.valDELTAE = str2double(VAP.settings.valDELTAE.Text);

%% Conditions
COND.valDENSITY = str2double(VAP.conditions.valDENSITY.Text);
VISC.valKINV = str2double(VAP.conditions.valKINV.Text);
try COND.valGUSTAMP = str2double(VAP.conditions.valGUSTAMP.Text); catch; COND.valGUSTAMP = 0; end
try COND.valGUSTL = str2double(VAP.conditions.valGUSTL.Text); catch; COND.valGUSTL = 0; end
try COND.valGUSTSTART = str2double(VAP.conditions.valGUSTSTART.Text); catch; COND.valGUSTSTART = 0; end

%% Vehicles
INPU.valVEHICLES = max(size(VAP.vehicle));

INPU.matVEHORIG = nan(INPU.valVEHICLES,3);
COND.vecVEHVINF = nan(INPU.valVEHICLES,1);
COND.vecVEHALPHA = nan(INPU.valVEHICLES,1);
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

VISC.vecFTURB = [];
VISC.vecFUSESECTIONS = [];
VISC.matFUSEAXIS = [];
VISC.matFUSEORIG = [];
VISC.matFGEOM = [];
VISC.matSECTIONFUSELAGE = [];
VISC.vecFUSEVEHICLE = [];
VISC.vecVSPANELS = [];
VISC.matVSGEOM = [];
VISC.vecFPANELS = [];
VISC.vecFPWIDTH = [];
VISC.vecINTERF = 0;

vecWINGINCID = [];
vecTRIMABLE = [];
vecWINGM = [];
matWINGORIG = [];
vecPANELS = [];
cellAIRFOILtemp = {};
VISC.cellAIRFOIL = {};
INPU.vecSYMtemp = [];
INPU.vecNtemp = [];
INPU.vecMtemp = [];
vecSECTIONS = [];
matSECTIONS = [];
vecSECTIONPANEL = [];

COND.vecWINGTRI = [];
COND.vecWAKETRI = [];

k = 1;
kk = 1;
kkk = 1;
kkkk = 1;

o = 1;
p = 1;

q = 1;

for i = 1:INPU.valVEHICLES
    
    try veh = VAP.vehicle{1,i}; catch; veh = VAP.vehicle; end
    
    INPU.matVEHORIG(i,:) = [str2double(veh.x.Text) str2double(veh.y.Text) str2double(veh.z.Text)];
    COND.vecVEHVINF(i,1) = str2double(veh.vinf.Text);
    try COND.vecVEHWEIGHT(i,1) = str2double(veh.weight.Text); catch; COND.vecVEHWEIGHT(i,1) = nan; end
    COND.vecVEHALPHA(i,1) = str2double(veh.alpha.Text);
    COND.vecVEHBETA(i,1) = str2double(veh.beta.Text);
    COND.vecVEHROLL(i,1) = str2double(veh.roll.Text);
    COND.vecVEHFPA(i,1) = str2double(veh.fpa.Text);
    COND.vecVEHTRK(i,1) = str2double(veh.trk.Text);
    
    try VEHI.vecVEHRADIUS(i,1) = str2double(veh.radius.Text); end
    
    try vecWINGS(i,1) = max(size(veh.wing)); catch; vecWINGS(i,1) = 0; end
    try vecSTRUCTURE(i,1) = max(size(veh.structure)); catch; vecSTRUCTURE(i,1) = 0; end
    try vecROTORS(i,1) = max(size(veh.rotor)); catch; vecROTORS(i,1) = 0; end
    try vecFUSELAGES(i,1) = max(size(veh.fuselage)); catch; vecFUSELAGES(i,1) = 0; end
    
    INPU.vecAREA(i) = str2double(veh.area.Text);
    INPU.vecSPAN(i) = str2double(veh.span.Text);
    INPU.vecCMAC(i) = str2double(veh.cmac.Text);
    
    %% Loading Wings
    for j = 1:vecWINGS(i)
        
        try win = veh.wing{1,j}; catch; win = veh.wing; end
        
        COND.vecWINGTRI(k,1) = nan;
        COND.vecWAKETRI(k,1) = nan;
        try if strcmpi(win.surfacetri.Text, 'true') COND.vecWINGTRI(k) = 1; end; end;
        try if strcmpi(win.waketri.Text, 'true') COND.vecWAKETRI(k) = 1; end; end;
        
        vecWINGINCID(k) = str2double(win.incidence.Text);
        if strcmpi(win.trimable.Text, 'true') vecTRIMABLE(j) = 1; else vecTRIMABLE(j) = 0; end
        
        vecWINGM(k,1) = str2double(win.M.Text);
        
        try matWINGORIG(k,:) = [str2double(win.xorig.Text) str2double(win.yorig.Text) str2double(win.zorig.Text)];
        catch; matWINGORIG(k,:) = [0 0 0]; end
        
        vecPANELS(k,1) = max(size(win.panel));
        
        for m = 1:vecPANELS(k,1)
            
            try pan = win.panel{1,m}; catch; pan = win.panel; end
            
            INPU.vecSYMtemp(kk,1) = floor(str2double(pan.symmetry.Text));
            try cellAIRFOILtemp{kk} = pan.airfoil.Text; end
            INPU.vecNtemp(kk,1) = floor(str2double(pan.N.Text));
            INPU.vecMtemp(kk,1) = floor(vecWINGM(k,1)); % Same for entire wing
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                sec = pan.section{1,n};
                
                matSECTIONS(kkk,:) = [str2double(sec.x.Text) + matWINGORIG(k,1) str2double(sec.y.Text) + matWINGORIG(k,2) str2double(sec.z.Text) + matWINGORIG(k,3) str2double(sec.chord.Text) vecWINGINCID(k)+str2double(sec.twist.Text) i];
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
    
    %% Loading Rotors
    for j = 1:vecROTORS(i)
        
        try rot = veh.rotor{1,j}; catch; rot = veh.rotor; end
        
        COND.vecROTORRPM(p,1) = str2double(rot.rpm.Text);
        INPU.vecROTDIAM(p,1) = str2double(rot.dia.Text);
        
        INPU.matROTORHUB(p,:) = [str2double(rot.xhub.Text) str2double(rot.yhub.Text) str2double(rot.zhub.Text)];
        INPU.matROTORAXIS(p,:) = [str2double(rot.axisx.Text) str2double(rot.axisy.Text) str2double(rot.axisz.Text)];
        
        INPU.vecROTORBLADES(p,:) = floor(str2double(rot.blades.Text));
        
        vecROTORM(p,1) = floor(str2double(rot.M.Text));
        
        try COND.vecCOLLECTIVE(p,1) = str2double(rot.collective.Text); catch; COND.vecCOLLECTIVE(p,1) = 0; end;
        
        vecPANELS(k,1) = max(size(rot.panel));
        
        flip = 1;
        try if strcmpi(rot.flipy.Text,'TRUE'); flip = -1; end; end
        
        for m = 1:vecPANELS(k,1)
            
            try pan = rot.panel{1,m}; catch; pan = rot.panel; end
            
            INPU.vecSYMtemp(kk,1) = 0;
            try cellAIRFOILtemp{kk} = pan.airfoil.Text; end
            INPU.vecNtemp(kk,1) = floor(str2double(pan.N.Text));
            INPU.vecMtemp(kk,1) = floor(vecROTORM(p,1)); % Same for entire wing
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                sec = pan.section{1,n};
                
                matSECTIONS(kkk,:) = [str2double(sec.x.Text) + INPU.matROTORHUB(p,1) flip*str2double(sec.y.Text) + INPU.matROTORHUB(p,2) str2double(sec.z.Text) + INPU.matROTORHUB(p,3) str2double(sec.chord.Text) str2double(sec.twist.Text) i];
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
    
    %% Fuselage
    
    for j = 1:vecFUSELAGES(i)
        try fuse = veh.fuselage{1,j}; catch; fuse = veh.fuselage; end
        
        VISC.vecFTURB(q,1) = str2double(fuse.transitionpanel.Text);
        
        VISC.vecFUSESECTIONS(q,1) = max(size(fuse.fsection));
        
        VISC.matFUSEAXIS(q,:) = [str2double(fuse.xaxis.Text) str2double(fuse.yaxis.Text) str2double(fuse.zaxis.Text)];
        VISC.matFUSEORIG(q,:) = [str2double(fuse.xorig.Text) str2double(fuse.yorig.Text) str2double(fuse.zorig.Text)];
        
        for m = 1:VISC.vecFUSESECTIONS(q,1)
            sec = fuse.fsection{1,m};
            
            VISC.matFGEOM(kkkk,:) = [str2double(sec.x.Text) str2double(sec.diameter.Text)];
            VISC.matSECTIONFUSELAGE(kkkk,1) = q;
            
            kkkk = kkkk + 1;
        end
        
        VISC.vecFUSEVEHICLE(q,1) = i;
        q = q + 1;
        
    end
end

INPU.valPANELS = sum(vecPANELS);


%% Reorganizing data for output

k = 1;

for i = 1:INPU.valPANELS
    
    sections = matSECTIONS(vecSECTIONPANEL == i,:);
    
    % (VAP3) add vehicle id to INPU.matGEOM
    %     sections(:,6) = VEHI.vecSURFACEVEHICLE(i);
    len = size(sections,1);
    
    if len == 2
        INPU.matGEOM(1,:,k) = sections(1,:);
        INPU.matGEOM(2,:,k) = sections(2,:);
        INPU.vecPANELWING(k,1) = INPU.vecPANELWING(i);
        INPU.vecPANELROTOR(k,1) = panel_rotors(i);
        
        if ~isempty(cellAIRFOILtemp)
            VISC.cellAIRFOIL{k,1} = cellAIRFOILtemp{i};
        end
        
        INPU.vecSYM(k,1) = INPU.vecSYMtemp(i);
        INPU.vecN(k,1) = INPU.vecNtemp(i);
        INPU.vecM(k,1) = INPU.vecMtemp(i);
        k = k + 1;
    else
        for j = 1:len - 1
            INPU.matGEOM(1,:,k) = sections(j,:);
            INPU.matGEOM(2,:,k) = sections(j+1,:);
            INPU.vecPANELWING(k,1) = INPU.vecPANELWING(i);
            INPU.vecPANELROTOR(k,1) = panel_rotors(i);
            
            INPU.vecN(k,1) = INPU.vecNtemp(i);
            INPU.vecM(k,1) = INPU.vecMtemp(i);
            
            if ~isempty(cellAIRFOILtemp)
                VISC.cellAIRFOIL{k,1} = cellAIRFOILtemp{i};
            end
            
            INPU.vecSYM(k,1) = 0;
            
            k = k + 1;
        end
        
        if INPU.vecSYMtemp(i) == 1
            INPU.vecSYM(k-len + 1) = 1;
        elseif INPU.vecSYMtemp(i) == 2
            INPU.vecSYM(k-1) = 2;
        end
        
    end
    
end

INPU.valPANELS = size(INPU.matGEOM,3);

if any(isnan(COND.vecVEHVINF))
    disp('One or more vehicle velocities was not read in - setting to unity and enabling fixed-lift analysis');
    COND.vecVEHVINF(isnan(COND.vecVEHVINF)) = ones(sum(isnan(COND.vecVEHVINF)),1);
    FLAG.FIXEDLIFT = 1;
end

%% Overwriting default values with user-specified
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








