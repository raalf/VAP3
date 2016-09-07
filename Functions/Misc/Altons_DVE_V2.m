% Jable Release
%   V0 - before fixing spanwise interp
%   V1 - fixed vertical panel (90deg dihedral)
%      - Reprogrammed the DVE interpolation method
%  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
%   V2 - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
%      - Comptue LE Sweep
%      - Project TE to DVE, Rotate adn Comptue TE Sweep

% INPUTS
clc
clear

showplot = 0;       %1 = show plot, 0 = hide plot
debug = 0;          %1 = clear all the variables at the end, Only show final matrices

% Panel Geometry
% Root Chord X, Y, Z, Chord, Twist
X1 = 0; Y1 = 0; Z1 = 0;   C1 = 2; e1 = 0;
X2 = 0; Y2 = 5; Z2 = 1; C2 = 1; e2 = 20;

% DVE Chordwise
n = 8; %chordwise

% DVE Spanwise
m = 3; %spanwise


% ==============================================================================

% Calculate twist at leading edge with dihedral correctly modelled
paneldata.rLE = [X1,Y1,Z1];
paneldata.tLE = [X2,Y2,Z2];
paneldata.rchord = C1;
paneldata.tchord = C2;

paneldata.repsilon = deg2rad(e1);
paneldata.tepsilon = deg2rad(e2);

% Read panel corners
panel = reshape(paneldesign(paneldata),3,4)';
panelX = [panel([1;4],1),panel([2;3],1)];
panelY = [panel([1;4],2),panel([2;3],2)];
panelZ = [panel([1;4],3),panel([2;3],3)];


% NChordwise and NSpanwise
% Generate extra points to find control point
N = n*2;
M = m*2;


% split Root Chord to chordwise elements
chordX(:,1) = linspace(panelX(1,1),panelX(2,1),M+1)';
chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
chordZ(:,1) = linspace(panelZ(1,1),panelZ(2,1),M+1)';
% split Tip Chord to chordwise elements
chordX(:,2) = linspace(panelX(1,2),panelX(2,2),M+1)';
chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
chordZ(:,2) = linspace(panelZ(1,2),panelZ(2,2),M+1)';

% linspace spanwise elemenets at each chordwise station

% Preallocate the memories for point matrices
PX2 = zeros(M+1,N+1);
PY2 = zeros(M+1,N+1);
PZ2 = zeros(M+1,N+1);

for i = 1:M+1
    if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
        spanwise = linspace(chordY(1,1),chordY(1,2),N+1);
        chordbase = chordY(1,:);
    else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
        spanwise = linspace(chordZ(1,1),chordZ(1,2),N+1);
        chordbase = chordZ(1,:);
    end
    PX2(i,:) = interp1(chordbase,chordX(i,:),spanwise); % X coordinates of all points
    PY2(i,:) = interp1(chordbase,chordY(i,:),spanwise); % Y coordinates of all points
    PZ2(i,:) = interp1(chordbase,chordZ(i,:),spanwise); % Z coordinates of all points
end



%% DVE Parameters Calculation
% Calculate Control Points, stored in 3D matrix
CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],m,n,3);
RightofCP = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],m,n,3);
LEmidpoint = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],m,n,3);
TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],m,n,3);
TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],m,n,3);


HS = elementwise3Ddist(RightofCP,CP);
HC = elementwise3Ddist(LEmidpoint,CP);
Area = HS.*HC.*4;

% Leading Edge Vectors
% Find leading edge vectors across entire panel
LE_vec = [PX2(1:2:end-1,2) PY2(1:2:end-1,2) PZ2(1:2:end-1,2)]-[PX2(1:2:end-1,1) PY2(1:2:end-1,1) PZ2(1:2:end-1,1)];
% Normalize LE vector
LE_norm = repmat(((LE_vec(:,1).^2+LE_vec(:,2).^2+LE_vec(:,3).^2).^0.5),1,3);
LE_vec = LE_vec./LE_norm;
% Repeat results n times for entire panel
LE_vec = repmat(reshape(LE_vec,m,1,3),1,n);


% Find Halfspan Vector. Vector from Control Point to LE midpoint
HS_vec = LEmidpoint - CP;
% Normalize Halfspan vector
HS_norm = (HS_vec(:,:,1).^2+HS_vec(:,:,2).^2+HS_vec(:,:,3).^2).^0.5;
HS_norm = repmat(HS_norm,1,1,3);
HS_vec = (HS_vec./HS_norm);

% Find Normal vector of the DVE.
% Leading Edge vector CROSS Halfspan vector
DVE_norm = cross(LE_vec,HS_vec,3);


% Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
roll = -atand(DVE_norm(:,:,2)./DVE_norm(:,:,3));

% Pitch in Degrees
% arcsin ( X component of DVE normal vector )
pitch  = asind(DVE_norm(:,:,1));

% Yaw = 0 Degrees
yaw = zeros(m,n);


% Rotate Reference Frame
% Leading Edge Vector in local reference frame
[ LE_local ] = glob_star_3D( LE_vec,roll,pitch,yaw );
% Compute DVE Leading Edge Sweep     arctan(LE X local component/ LE Y local component)
Sweep_LE = atand(LE_local(:,:,1)./LE_local(:,:,2));

% % Compute DVE Trailing Edge Sweep
% Project TE Points onto DVE plane 
% (TE_Left / TE_Right) (CP)                   (DVE_norm)
% q(x,y,z) TE point | p(a,b,c) Control Point | n(d,e,f) DVE normal
% q_proj = q - dot(q-p,n)*n
TE_Left_proj = TE_Left-repmat(dot(TE_Left-CP,DVE_norm,3),1,1,3).*DVE_norm;
TE_Right_proj = TE_Right-repmat(dot(TE_Right-CP,DVE_norm,3),1,1,3).*DVE_norm;
TE_vec_proj = TE_Right_proj - TE_Left_proj;

% Rotate the Projected TE on DVE to local reference frame
[ TE_local ] = glob_star_3D( TE_vec_proj,roll,pitch,yaw );
Sweep_TE = atand(TE_local(:,:,1)./TE_local(:,:,2));


% Compute DVE Mid-chord Sweep
% Average of LE and TE Sweep
Sweep_MID = (Sweep_LE+Sweep_TE)./2;

count = length(Area(:));

DVE.xo = reshape(CP(:),count,3);
DVE.xsi = HC(:);
DVE.eta = HS(:);
DVE.phiLE = Sweep_LE(:);
DVE.phiTE = Sweep_TE(:);
DVE.roll = roll(:);
DVE.pitch = pitch(:);
DVE.yaw = yaw(:);
DVE.Index = 1:count;

break

%%
%ALL PLOTS
if showplot == 1
    clf
    figure(1)
    hold on
    mesh(PX2(1:2:end,1:2:end),PY2(1:2:end,1:2:end),PZ2(1:2:end,1:2:end));
    %     colormap white
    scatter3(reshape(CP(:,:,1),n*m,1),reshape(CP(:,:,2),n*m,1),reshape(CP(:,:,3),n*m,1),'xr');
    
    % TE Points
    scatter3(reshape(TE_Right(:,:,1),n*m,1),reshape(TE_Right(:,:,2),n*m,1),reshape(TE_Right(:,:,3),n*m,1),'rx');
    scatter3(reshape(TE_Left(:,:,1),n*m,1),reshape(TE_Left(:,:,2),n*m,1),reshape(TE_Left(:,:,3),n*m,1),'go');
    scatter3(reshape(TE_Left_proj(:,:,1),n*m,1),reshape(TE_Left_proj(:,:,2),n*m,1),reshape(TE_Left_proj(:,:,3),n*m,1),'ks');

    
%     
%     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),HS_vec(:,:,1),HS_vec(:,:,2),HS_vec(:,:,3));
%     quiver3(LEmidpoint(:,:,1),LEmidpoint(:,:,2),LEmidpoint(:,:,3),LE_vec(:,:,1),LE_vec(:,:,2),LE_vec(:,:,3));
%     % DVE Normal Vector
%     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),DVE_norm(:,:,1),DVE_norm(:,:,2),DVE_norm(:,:,3));
%     
%     quiver3(LEmidpoint(:,:,1),LEmidpoint(:,:,2),LEmidpoint(:,:,3),LE_local(:,:,1),LE_local(:,:,2),LE_local(:,:,3));
%     
    quiver3(TE_Left_proj(:,:,1),TE_Left_proj(:,:,2),TE_Left_proj(:,:,3),TE_vec_proj(:,:,1),TE_vec_proj(:,:,2),TE_vec_proj(:,:,3));

    %     plot3(panel(:,1),panel(:,2),panel(:,3))
    %     surf(panelX,panelY,panelZ)
    hold off
    axis equal
    grid on
end

if debug~=1
    clear chordX chordY chordZ panelX panelY panelZ debug
    clear LE_norm HC HS
    clear X1 X2 Y1 Y2 Z1 Z2 AX1 AZ1 AX2 AZ2 C1 C2 e2 e1 chordbase
    clear panel paneldata i j M N
    clear  halfchord halfspan
end
