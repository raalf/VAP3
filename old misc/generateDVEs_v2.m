function [Aircraft, FW] = generateDVEs_v2(Aircraft, FW, ii)
%   V0 - before fixing spanwise interp
%   V1 - fixed vertical panel (90deg dihedral)
%      - Reprogrammed the DVE interpolation method
%  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
%   V2 - Rework the Leading Edge Vector
%      - Caculate the Yaw angle of the DVE by assuming yaw=0 to rotate the xsi vector
%      - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
%      - Comptue LE Sweep
%      - Project TE to DVE, Rotate adn Comptue TE Sweep

% Fixed how DVEs matrix is converted from 2D grid to 1D array. 16/01/2016 (Alton)

showplot = 0;       %1 = show plot, 0 = hide plot
debug = 0;          %1 = clear all the variables at the end, Only show final matrices

Panel = Aircraft.Surface(ii).Panel;

X1 = Panel.X(1); Y1 = Panel.Y(1); Z1 = Panel.Z(1); C1 = Panel.Chord(1); e1 = Panel.Twist(1);
X2 = Panel.X(2); Y2 = Panel.Y(2); Z2 = Panel.Z(2); C2 = Panel.Chord(2); e2 = Panel.Twist(2);

n = FW.Panels(ii).n;
% m = FW.Panels(ii).m;
m = FW.m;

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

panelTE = [panel(end,:); panel(end-1,:)];


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
CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],m,n,3);
LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],m,n,3);
TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],m,n,3);
TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],m,n,3);

% Filter Leading Edge Point
LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],m,n,3);
LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],m,n,3);

% Create eta vector for full leading edge
% Non-normalized
LE_vec = LE_Right - LE_Left;

% Create half chord xsi vector
% Non-normalized
xsi_vec = LE_Mid - CP;


%
DVE_norm = normalize3D(cross(LE_vec, xsi_vec, 3));



% Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
% atan2d is used here
% roll(nu) right wing up positive
nu = -atan2d(DVE_norm(:,:,2),DVE_norm(:,:,3));

% Pitch in Degrees
% arcsin ( X component of DVE normal vector )
epsilon = asind(DVE_norm(:,:,1));

% Yaw in Degrees
% xsi in local with roll picth, yaw set to zero.. but WHY?
xsi_local = glob_star_3D( xsi_vec,nu,epsilon,zeros(m,n) );
% Magnitude of half chord vector
xsi = (xsi_local(:,:,1).^2+xsi_local(:,:,2).^2+xsi_local(:,:,3).^2).^0.5;
psi = atand(xsi_local(:,:,2)./xsi_local(:,:,1));

% Find eta. bring non-normalized LE_vec to local and half the Y component
LE_vec_local = glob_star_3D( LE_vec,nu,epsilon,psi);
eta = LE_vec_local(:,:,2)./2;



% Find Leading Edge Sweep
% arctan(LE X local component/ LE Y local component)
phi_LE = atand(LE_vec_local(:,:,1)./LE_vec_local(:,:,2));



% Find Trailing Edge Sweep
    % Project TE Points onto DVE plane
    % (TE_Left / TE_Right) (CP)                   (DVE_norm)
    % q(x,y,z) TE point | p(a,b,c) Control Point | n(d,e,f) DVE normal
    % q_proj = q - dot(q-p,n)*n
    TE_Left_proj = TE_Left-repmat(dot(TE_Left-CP,DVE_norm,3),1,1,3).*DVE_norm;
    TE_Right_proj = TE_Right-repmat(dot(TE_Right-CP,DVE_norm,3),1,1,3).*DVE_norm;
    TE_vec_proj = TE_Right_proj - TE_Left_proj;

% Rotate the Projected TE on DVE to local reference frame
% arctan(Projected TE local X component/Projected TE local Y component)
TE_vec_proj_local = glob_star_3D( TE_vec_proj,nu,epsilon,psi );
phi_TE = atand(TE_vec_proj_local(:,:,1)./TE_vec_proj_local(:,:,2));


% Compute DVE Mid-chord Sweep
% Average of LE and TE Sweep
phi_MID = (phi_LE+phi_TE)./2;

% Calculating Area
Area = eta.*xsi.*4;

% Calculating local free stream velocity at xo location


%% Formatting the output

% out = [Name eta xsi nu eps psi xo(x y z) phiLE phiMID phiTE Area norm(x y z) u(x y z)]

name = [1:(n*m)]';

if ii > 1
    name = name + sum([FW.Panels(1:ii-1).n])*m;
end

FW.Panels(ii).DVE.Index = name;
% FW.Panels(ii).DVE.area = reshape(Area',[],1);
% FW.Panels(ii).DVE.eta = reshape(eta,[],1);
% FW.Panels(ii).DVE.xsi = reshape(xsi,[],1);
% FW.Panels(ii).DVE.roll = reshape(nu,[],1);
% FW.Panels(ii).DVE.pitch = reshape(epsilon,[],1);
% FW.Panels(ii).DVE.yaw = reshape(psi,[],1);
% FW.Panels(ii).DVE.xo = [reshape(CP(:,:,1)',[],1) reshape(CP(:,:,2)',[],1) reshape(CP(:,:,3)',[],1)];
% FW.Panels(ii).DVE.phiLE = reshape(phi_LE',[],1);
% FW.Panels(ii).DVE.phiMID = reshape(phi_MID',[],1);
% FW.Panels(ii).DVE.phiTE = reshape(phi_TE',[],1);
% FW.Panels(ii).DVE.norm = [reshape(DVE_norm(:,:,1)',[],1) reshape(DVE_norm(:,:,2)',[],1) reshape(DVE_norm(:,:,3)',[],1)];
FW.Panels(ii).DVE.singfct = 0;
% FW.Panels(ii).DVE.TECoordL = [reshape(TE_Left_proj(:,:,1)',[],1) reshape(TE_Left_proj(:,:,2)',[],1) reshape(TE_Left_proj(:,:,3)',[],1)];
% FW.Panels(ii).DVE.TECoordR = [reshape(TE_Right_proj(:,:,1)',[],1) reshape(TE_Right_proj(:,:,2)',[],1) reshape(TE_Right_proj(:,:,3)',[],1)];
FW.Panels(ii).TECoord = panelTE;

%% Version2 -> Reformat 2D vector to 1D
count = length(Area(:));

FW.Panels(ii).DVE.area = reshape(Area',count,1);%Area(:);
FW.Panels(ii).DVE.eta = reshape(eta',count,1);%eta(:);
FW.Panels(ii).DVE.xsi = reshape(xsi',count,1);%xsi(:);
FW.Panels(ii).DVE.roll = reshape(nu',count,1);%nu(:);
FW.Panels(ii).DVE.pitch = reshape(epsilon',count,1);%epsilon(:);
FW.Panels(ii).DVE.yaw = reshape(psi',count,1);%psi(:);
FW.Panels(ii).DVE.xo = reshape(permute(CP, [2 1 3]),count,3);%reshape(CP(:),count,3);
FW.Panels(ii).DVE.phiLE = reshape(phi_LE',count,1);%phi_LE(:);
FW.Panels(ii).DVE.phiMID = reshape(phi_MID',count,1);%phi_MID(:);
FW.Panels(ii).DVE.phiTE = reshape(phi_TE',count,1);%phi_TE(:);
FW.Panels(ii).DVE.norm = reshape(permute(DVE_norm, [2 1 3]),count,3);%reshape(DVE_norm(:),count,3);
FW.Panels(ii).DVE.TECoordL = reshape(permute(TE_Left_proj, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
FW.Panels(ii).DVE.TECoordR = reshape(permute(TE_Right_proj, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);
FW.Panels(ii).DVE.LECoordL = reshape(permute(LE_Left, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
FW.Panels(ii).DVE.LECoordR = reshape(permute(LE_Right, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);

%%

FW.Panels(ii).Edge1 = [FW.Panels(ii).DVE.Index(1):FW.Panels(ii).n:FW.Panels(ii).DVE.Index(end)]'; % Getting all elements in the chordwise direction (1 x m vector) at EDGE 1. This will come in handy later on when setting BC between panels
FW.Panels(ii).Edge2 = FW.Panels(ii).Edge1+FW.Panels(ii).n-1; % Getting all elements in the chordwise direction at EDGE 2 by adding the number of elements in the spanwise direction - 1



%% ALL PLOTS
if showplot == 1
    clf
    figure(1)
    hold on
    axis equal
    surf(PX2(1:2:end,1:2:end),PY2(1:2:end,1:2:end),PZ2(1:2:end,1:2:end));
    colormap white
    scatter3(reshape(CP(:,:,1),n*m,1),reshape(CP(:,:,2),n*m,1),reshape(CP(:,:,3),n*m,1),'xg');
    
%     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),HS_vec(:,:,1),HS_vec(:,:,2),HS_vec(:,:,3));
%     quiver3(LEmidpoint(:,:,1),LEmidpoint(:,:,2),LEmidpoint(:,:,3),LE_vec(:,:,1),LE_vec(:,:,2),LE_vec(:,:,3));
    quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),DVE_norm(:,:,1),DVE_norm(:,:,2),DVE_norm(:,:,3));
    
    % plot3(panel(:,1),panel(:,2),panel(:,3))
    % surf(panelX,panelY,panelZ)
    hold off
end

%% Debug Mode
if debug~=1
    clear chordX chordY chordZ panelX panelY panelZ m n debug
    clear X1 X2 Y1 Y2 Z1 Z2 AX1 AZ1 AX2 AZ2 C1 C2 e2
    clear panel paneldata i j M N
    clear  halfchord halfspan
end

end















