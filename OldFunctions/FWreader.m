function [FW, Aircraft, Conditions] = FWreader(file)

% clc
% clear

fp = fopen(file);

%% Reading header flags
% Reading relaxed wake flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Relax = fscanf(fp,'%d');

% Reading steady or unsteady flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Steady = fscanf(fp,'%d');

% Reading symmetry flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Sym = fscanf(fp,'%d');

% Reading trim flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Trim = fscanf(fp,'%d');
%% Reading time step information
% Reading maximum number of time steps
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Maxtime = fscanf(fp,'%d');

% Reading time step width
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Deltime = fscanf(fp,'%lf');

% Reading deltae
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Deltae = fscanf(fp,'%lf');
%% Reading flow conditions
% Reading freestream velocity
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Uinf = fscanf(fp,'%lf');

% Reading alpha start, stop and increment into column vector
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Conditions.Alpha = fscanf(fp,'%lf');

% Reading sideslip angle
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Conditions.Beta = fscanf(fp,'%lf');

% Reading density
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Conditions.Density = fscanf(fp,'%lf');

% Reading density
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Conditions.KinematicViscosity = fscanf(fp,'%lf');
%% Reading Aircraft Reference Values
% Reading wing area
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.Reference.S = fscanf(fp,'%lf');

% Reading wing span
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.Reference.b = fscanf(fp,'%lf');

% Reading mean aerodynamic chord
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.Reference.Cmac = fscanf(fp,'%lf');

% Reading aircraft weight
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.Reference.Weight = fscanf(fp,'%lf');

% Reading CG location into column vector
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.Reference.CG = fscanf(fp,'%lf');

% Reading CMo of wing
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.General.CMo = fscanf(fp,'%lf');
%% Reading panel/wing/lifting line information
% Reading No. of Wings
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.General.Wings = fscanf(fp,'%lf');

% Reading No. of panels
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
Aircraft.General.Panels = fscanf(fp,'%lf');

% Reading No. of chordwise lifting lines
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.m = fscanf(fp,'%lf');

% Reading No. of airfoils
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Airfoils = fscanf(fp,'%lf');
%% Reading panel information and geometry

FW.n = zeros(Aircraft.General.Panels,1);
Aircraft.General.Airfoil = zeros(Aircraft.General.Panels,1);

for i = 1:Aircraft.General.Panels
    % Reading 'n'
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    FW.n(i) = fscanf(fp,'%lf',1);
    
    % Reading airfoil number
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    Aircraft.General.Airfoil(i) = fscanf(fp,'%lf',1);
    
    % Reading left panel number
    ch = fscanf(fp,'%c',1);
    while(ch~=':');
        ch = fscanf(fp,'%c',1);
    end
    FW.Left(i) = fscanf(fp,'%lf',1);
    
    % Reading left panel number
    ch = fscanf(fp,'%c',1);
    while(ch~=':');
        ch = fscanf(fp,'%c',1);
    end
    FW.Right(i) = fscanf(fp,'%lf',1);
    
    % Skipping geometry column header
    fgets(fp);
    fgets(fp);
    
    % Reading geometry
    % Explanation below:
    %{
        info_geometry(x,y,z)
            x is for the left or right point
                1 left
                2 right
            y is for the values
                1 x
                2 y
                3 z
                4 chord
                5 epsilon
                6 boundary condition
            z is panel number
    %}
    
    FW.Geometry(1,:,i) = fscanf(fp,'%lf');
    fgets(fp); % Skipping another column header
    FW.Geometry(2,:,i) = fscanf(fp,'%lf');
    
end

FW.BC = reshape(FW.Geometry(:,end,:),[],1);

%% Reading vertical tail information
% Reading number of panels
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.VerticalStabilizer.Panels = fscanf(fp,'%lf',1);

% Skipping vtail geometry column headers
fgets(fp);
fgets(fp);

FW.VerticalStabilizer.Geometry = zeros(FW.VerticalStabilizer.Panels,4);

% Reading v-stab geometry
% Explanation below:
%{
    info_vgeometry(x,y)
        x is the panel number
        y is for the values
            1 panel number
            2 panel chord
            3 panel area
            4 panel airfoil
%}

for j = 1:FW.VerticalStabilizer.Panels
    FW.VerticalStabilizer.Geometry(j,:) = fscanf(fp,'%lf');
end
%% Reading fuselage information
% Reading number of sections
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Fuselage.Sections = fscanf(fp,'%lf',1);

% Reading width of sections
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Fuselage.Width = fscanf(fp,'%lf',1);

% Reading turbulence transition point
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
FW.Fuselage.Turbulent = fscanf(fp,'%lf',1);

% Skipping fuselage geometry column headers (I don't know why I need 3)
fgetl(fp);
fgetl(fp);
% fgets(fp);

FW.Fuselage.Geometry = zeros(FW.Fuselage.Sections,2);

for j = 1:FW.Fuselage.Sections
    FW.Fuselage.Geometry(j,:) = fscanf(fp,'%lf',2);
    % I have no idea why I need these:
    %     fgetl(fp);
    %     fgetl(fp);
    % Without them, fscanf returns an extra section number which isn't in
    % the text file
end

% Reading intereference drag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end

FW.Interference = fscanf(fp,'%lf',1);

fclose(fp);

%% Putting Panels into Structures
for ii = 1:Aircraft.General.Panels
    
    Aircraft.Surface(ii).Name = sprintf('Panel %d',ii);
    Aircraft.Surface(ii).X = 0;
    Aircraft.Surface(ii).Y = 0;
    Aircraft.Surface(ii).Z = 0;
    Aircraft.Surface(ii).IncidenceAngle = 0;
   
    
    FW.Panels(ii).Name = sprintf('Panel %d',ii);
    FW.Panels(ii).n = FW.n(ii);
    FW.Panels(ii).m = FW.m;
    FW.Panels(ii).Airfoil = Aircraft.General.Airfoil(ii);
    
    Aircraft.Surface(ii).Panel.X(1) = FW.Geometry(1,1,ii);
    Aircraft.Surface(ii).Panel.Y(1) = FW.Geometry(1,2,ii);
    Aircraft.Surface(ii).Panel.Z(1) = FW.Geometry(1,3,ii);
    Aircraft.Surface(ii).Panel.Chord(1) = FW.Geometry(1,4,ii);
    Aircraft.Surface(ii).Panel.Twist(1) = FW.Geometry(1,5,ii);
    Aircraft.Surface(ii).Panel.X(2) = FW.Geometry(2,1,ii);
    Aircraft.Surface(ii).Panel.Y(2) = FW.Geometry(2,2,ii);
    Aircraft.Surface(ii).Panel.Z(2) = FW.Geometry(2,3,ii);
    Aircraft.Surface(ii).Panel.Chord(2) = FW.Geometry(2,4,ii);
    Aircraft.Surface(ii).Panel.Twist(2) = FW.Geometry(2,5,ii);
    
end

% These are used for determining the number of spanwise elements in each
% wing

Aircraft.Reference.AR = Aircraft.Reference.b^2/Aircraft.Reference.S; % Reference aspect ratio

FW = rmfield(FW, {'BC','Left','Right','Trim','Airfoils','n','Geometry'});
Aircraft.General = rmfield(Aircraft.General, {'Airfoil'});

clear ans ch fp i ii j temp wing_start_index


end




















