function [...
    flagRELAX, flagSTEADY, flagTRI, matGEOM, valMAXTIME, valMINTIME, valDELTIME, valDELTAE, ...
    valDENSITY, valKINV, valVEHICLES, matVEHORIG, vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHROLL, ...
    vecVEHFPA, vecVEHTRK, vecWINGS, vecWINGINCID, vecAREA, vecSPAN, vecCMAC, vecWINGM, ...
    vecPANELS, vecSYM, vecN, vecM, vecSECTIONS, matSECTIONS, vecSECTIONPANEL, vecWING, ...
    vecWINGVEHICLE, valPANELS, vecROTORS, vecROTORRPM, vecROTDIAM, matROTORHUB, vecROTORAXIS, vecROTORBLADES,...
    vecROTORM, vecROTOR, vecFTURB, vecFUSESECTIONS, matFGEOM, matSECTIONFUSELAGE, vecFUSEVEHICLE, matFUSEAXIS, matFUSEORIG...
    ] = fcnXMLREAD(filename)

% clc
% clear
%
% filename = 'inputs/XMLtest.vap';

% OUTPUT

% vecROTORS - Rows are wing number, tells which vehicle has how many rotors
% vecROTORRPM - Rows are rotor number, tells the rotor RPM
% vecROTDIAM - Rows are rotor number, tells the rotor diameter
% vecROTORHUB - Rows are rotor number, column is xyz in vehicle local. Tells rotor hub location
% vecROTORAXIS - Rows are rotor number, column is xyz in vehicle local. Unit vector of rotor plane normal
% vecROTORM - Rows are rotor number, tells how many chordwise elements on the rotor blade
% vecROTOR - Rows are panel number, tells which rotor the panel belongs to

% matSURFACETYPE - surface x 2 matrix defining which surface is a wing or rotor, where 1 in column 1 is wing, 1 in column 2 is rotor

% vecWINGVEHICLE - Rows are wing number, tells which vehicle each wing belongs to
% vecWING - Rows are panel number, tells which wing each panel belongs to
% valPANELS - Total number of panels
% vecPANELS - Rows are wing number, tells how many panels are on each wing
% vecWINGM - Rows are wing number, tells us the number of chordwise lifting lines on each wing
% vecWINGCMAC - Rows are wing number, tells us the mean aerodynamic chord of each wing
% vecWINGSPAN - Rows are wing number, tell us the span of each wing (tip to tip)
% vecSYM - Rows are panel number, tells us whether edge 1 or 2 of the panel has symmetry boundary condition
% vecN - Rows are panel number, tells us the number of spanwise DVEs for each panel
% vecM - Rows are panel number, tells us the number of chordwise DVEs for each panel (same as vecWINGM but on the DVE level)
% vecWINGINCID - Rows are wing number, tells us the incidence angle of the wing
% matSECTIONS - number of sections x 5 matrix of [x y z chord twist] for each section

%%
inp = fcnXML2STRUCT(filename);
VAP = inp.VAP;

%% Settings
if strcmpi(VAP.settings.flagRELAX.Text, 'true') flagRELAX = 1; else flagRELAX = 0; end
if strcmpi(VAP.settings.flagSTEADY.Text, 'true') flagSTEADY = 1; else flagSTEADY = 0; end
if strcmpi(VAP.settings.flagTRI.Text, 'true') flagTRI = 1; else flagTRI = 0; end

valMAXTIME = floor(str2double(VAP.settings.valMAXTIME.Text));
valMINTIME = floor(str2double(VAP.settings.valMINTIME.Text));
valDELTIME = str2double(VAP.settings.valDELTIME.Text);
valDELTAE = str2double(VAP.settings.valDELTAE.Text);

%% Conditions
valDENSITY = str2double(VAP.conditions.valDENSITY.Text);
valKINV = str2double(VAP.conditions.valKINV.Text);

%% Vehicles
valVEHICLES = max(size(VAP.vehicle));

matVEHORIG = nan(valVEHICLES,3);
vecVEHVINF = nan(valVEHICLES,1);
vecVEHALPHA = nan(valVEHICLES,1);
vecVEHBETA = nan(valVEHICLES,1);
vecVEHROLL = nan(valVEHICLES,1);
vecVEHFPA = nan(valVEHICLES,1);
vecVEHTRK = nan(valVEHICLES,1);

vecWINGS = nan(valVEHICLES,1);
vecROTORS = nan(valVEHICLES,1);
vecFUSELAGES = nan(valVEHICLES,1);

vecROTORRPM = [];
vecROTDIAM = [];
matROTORHUB = [];
vecROTORAXIS = [];
vecROTORBLADES = [];
vecROTORM = [];

vecFTURB = [];
vecFUSESECTIONS = [];
matFUSEAXIS = [];
matFUSEORIG = [];
matFGEOM = [];
matSECTIONFUSELAGE = [];
vecFUSEVEHICLE = [];

vecWINGINCID = [];
vecTRIMABLE = [];
vecWINGM = [];
matWINGORIG = [];
vecPANELS = [];
vecSYMtemp = [];
vecNtemp = [];
vecMtemp = [];
vecSECTIONS = [];
matSECTIONS = [];
vecSECTIONPANEL = [];

k = 1;
kk = 1;
kkk = 1;
kkkk = 1;

p = 1;

q = 1;

for i = 1:valVEHICLES
    
    try veh = VAP.vehicle{1,i}; catch; veh = VAP.vehicle; end
    
    matVEHORIG(i,:) = [str2double(veh.x.Text) str2double(veh.y.Text) str2double(veh.z.Text)];
    vecVEHVINF(i,1) = str2double(veh.vinf.Text);
    vecVEHALPHA(i,1) = str2double(veh.alpha.Text);
    vecVEHBETA(i,1) = str2double(veh.beta.Text);
    vecVEHROLL(i,1) = str2double(veh.roll.Text);
    vecVEHFPA(i,1) = str2double(veh.fpa.Text);
    vecVEHTRK(i,1) = str2double(veh.trk.Text);
    
    try vecWINGS(i,1) = max(size(veh.wing)); catch; vecWINGS(i,1) = 0; end
    try vecROTORS(i,1) = max(size(veh.rotor)); catch; vecROTORS(i,1) = 0; end
    try vecFUSELAGES(i,1) = max(size(veh.fuselage)); catch; vecFUSELAGES(i,1) = 0; end
    
    vecAREA(i) = str2double(veh.area.Text);
    vecSPAN(i) = str2double(veh.span.Text);
    vecCMAC(i) = str2double(veh.cmac.Text);
    
    %% Loading Wings
    for j = 1:vecWINGS(i)
              
        try win = veh.wing{1,j}; catch; win = veh.wing; end
        
        vecWINGINCID(k) = str2double(win.incidence.Text);
        if strcmpi(win.trimable.Text, 'true') vecTRIMABLE(j) = 1; else vecTRIMABLE(j) = 0; end
        
        vecWINGM(k,1) = str2double(win.M.Text);
        
        try matWINGORIG(k,:) = [str2double(win.xwingorig.Text) str2double(win.ywingorig.Text) str2double(win.zwingorig.Text)];
        catch; matWINGORIG(k,:) = [0 0 0]; end
        
        vecPANELS(k,1) = max(size(win.panel));
        
        for m = 1:vecPANELS(k,1)
            
            try pan = win.panel{1,m}; catch; pan = win.panel; end
            
            vecSYMtemp(kk,1) = floor(str2double(pan.symmetry.Text));
            vecNtemp(kk,1) = floor(str2double(pan.N.Text));
            vecMtemp(kk,1) = floor(vecWINGM(k,1)); % Same for entire wing
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                sec = pan.section{1,n};
                
                matSECTIONS(kkk,:) = [str2double(sec.x.Text) + matWINGORIG(k,1) str2double(sec.y.Text) + matWINGORIG(k,2) str2double(sec.z.Text) + matWINGORIG(k,3) str2double(sec.chord.Text) vecWINGINCID(k)+str2double(sec.twist.Text) i];
                vecSECTIONPANEL(kkk,1) = kk;
                
                kkk = kkk + 1;
            end
            
            vecPANELWING(kk,1) = k;
            vecPANELROTOR(kk,1) = 0;
            
            kk = kk + 1;
        end
        
        vecWINGVEHICLE(k,1) = i;
        k = k + 1;
    end
    
    %% Loading Rotors
    for j = 1:vecROTORS(i)
        
        try rot = veh.rotor{1,j}; catch; rot = veh.rotor; end
        
        vecROTORRPM(p,1) = str2double(rot.rpm.Text);
        vecROTDIAM(p,1) = str2double(rot.dia.Text);
        
        matROTORHUB(p,:) = [str2double(rot.xhub.Text) str2double(rot.yhub.Text) str2double(rot.zhub.Text)];
        vecROTORAXIS(p,:) = [str2double(rot.axisx.Text) str2double(rot.axisy.Text) str2double(rot.axisz.Text)];
        
        vecROTORBLADES(p,:) = floor(str2double(rot.blades.Text));
        
        vecROTORM(p,1) = floor(str2double(rot.M.Text));
        
        vecPANELS(k,1) = max(size(rot.panel));
        
        flip = 1;
        try if strcmpi(rot.flipy.Text,'TRUE'); flip = -1; end; end
        
        for m = 1:vecPANELS(k,1)
            
            try pan = rot.panel{1,m}; catch; pan = rot.panel; end
            
            vecSYMtemp(kk,1) = 0;
            vecNtemp(kk,1) = floor(str2double(pan.N.Text));
            vecMtemp(kk,1) = floor(vecROTORM(p,1)); % Same for entire wing
            
            vecSECTIONS(kk,1) = max(size(pan.section));
            
            for n = 1:vecSECTIONS(kk,1)
                sec = pan.section{1,n};
                
                matSECTIONS(kkk,:) = [str2double(sec.x.Text) + matROTORHUB(p,1) flip*str2double(sec.y.Text) + matROTORHUB(p,2) str2double(sec.z.Text) + matROTORHUB(p,3) str2double(sec.chord.Text) str2double(sec.twist.Text) i];
                vecSECTIONPANEL(kkk,1) = kk;
                
                kkk = kkk + 1;
            end
            
            vecPANELWING(kk,1) = 0;
            vecPANELROTOR(kk,1) = p;
            
            kk = kk + 1;
        end
        
        p = p + 1;
        
        vecWINGVEHICLE(k,1) = i;
        k = k + 1;
    end
    
    %% Fuselage
   
    for j = 1:vecFUSELAGES(i)
        try fuse = veh.fuselage{1,j}; catch; fuse = veh.fuselage; end
        
        vecFTURB(q,1) = str2double(fuse.transitionpanel.Text);
        
        vecFUSESECTIONS(q,1) = max(size(fuse.fsection));
        
        matFUSEAXIS(q,:) = [str2double(fuse.xaxis.Text) str2double(fuse.yaxis.Text) str2double(fuse.zaxis.Text)];
        matFUSEORIG(q,:) = [str2double(fuse.xorig.Text) str2double(fuse.yorig.Text) str2double(fuse.zorig.Text)];
        
        for m = 1:vecFUSESECTIONS(q,1)
            sec = fuse.fsection{1,m};
            
            matFGEOM(kkkk,:) = [str2double(sec.x.Text) str2double(sec.diameter.Text)];
            matSECTIONFUSELAGE(kkkk,1) = q;
            
            kkkk = kkkk + 1; 
        end

        vecFUSEVEHICLE(q,1) = i;
        q = q + 1;
        
    end
end

valPANELS = sum(vecPANELS);


%% Reorganizing data for output

k = 1;

for i = 1:valPANELS
    
    sections = matSECTIONS(vecSECTIONPANEL == i,:);
    
    % (VAP3) add vehicle id to matGEOM
%     sections(:,6) = vecWINGVEHICLE(i);
    len = size(sections,1);
    
    if len == 2
        matGEOM(1,:,k) = sections(1,:);
        matGEOM(2,:,k) = sections(2,:);
        vecWING(k,1) = vecPANELWING(i);
        vecROTOR(k,1) = vecPANELROTOR(i);
        vecSYM(k,1) = vecSYMtemp(i);
        vecN(k,1) = vecNtemp(i);
        vecM(k,1) = vecMtemp(i);
        k = k + 1;
    else
        for j = 1:len - 1
            matGEOM(1,:,k) = sections(j,:);
            matGEOM(2,:,k) = sections(j+1,:);
            vecWING(k,1) = vecPANELWING(i);
            vecROTOR(k,1) = vecPANELROTOR(i);
            
            vecN(k,1) = vecNtemp(i);
            vecM(k,1) = vecMtemp(i);
            vecSYM(k,1) = 0;
            
            k = k + 1;
        end
        
        if vecSYMtemp(i) == 1
            vecSYM(k-len) = 1;
        elseif vecSYMtemp(i) == 2
            vecSYM(k-1) = 2;
        end
        
    end
    
end



%% Duplicating rotor blades
% THIS SHIT DON'T WORK
% [matNEWGEOM, ~, vecNnew, vecMnew, vecSYMnew] = fcnMULTIBLADE(sum(vecROTOR > 0), vecROTORBLADES, vecN(vecROTOR > 0), vecM(vecROTOR > 0), vecSYM(vecROTOR > 0), matGEOM(:,:,vecROTOR > 0));
%
% matGEOM = cat(3, matGEOM, matNEWGEOM);
%
% vecN = cat(1, vecN, vecNnew);
% vecM = cat(1, vecM, vecMnew);
% vecSYM = cat(1, vecSYM, vecSYMnew);

valPANELS = size(matGEOM,3);



end








