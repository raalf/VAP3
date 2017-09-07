% NOT Jable Release
%   V0 - before fixing spanwise interp
%   V1 - fixed vertical panel (90deg dihedral)
%      - Reprogrammed the DVE interpolation method
%  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
%   V2 - Rework the Leading Edge Vector
%      - Caculate the Yaw angle of the DVE by assuming yaw=0 to rotate the xsi vector
%      - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
%      - Comptue LE Sweep
%      - Project TE to DVE, Rotate adn Comptue TE Sweep

% INPUTS
clc
clear

showplot = 1;       %1 = show plot, 0 = hide plot
debug = 0;          %1 = clear all the variables at the end, Only show final matrices

% Panel Geometry
% Root Chord X, Y, Z, Chord, Twist
X1 = -1.349; Y1 = 0.09069; Z1 = -0.1778; e1 = -20; C1 = 2.0;
X2 = 8.956  ;Y2 = 1.959; Z2 = 1.953;   e2 = 50; C2 = 0.4;

% DVE Size
n = 10; %spanwise
m = 5; %chordwise


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







%%


%ALL PLOTS
if showplot == 1
    clf
    figure(1)
    hold on
%     mesh(PX2(1:2:end,1:2:end),PY2(1:2:end,1:2:end),PZ2(1:2:end,1:2:end));
    %     colormap white
%     scatter3(reshape(CP(:,:,1),n*m,1),reshape(CP(:,:,2),n*m,1),reshape(CP(:,:,3),n*m,1),'xr');
    
    % TE Points
    
    
%     scatter3(reshape(LE_Left(:,:,1),n*m,1),reshape(LE_Left(:,:,2),n*m,1),reshape(LE_Left(:,:,3),n*m,1),'b');
%     scatter3(reshape(LE_Right(:,:,1),n*m,1),reshape(LE_Right(:,:,2),n*m,1),reshape(LE_Right(:,:,3),n*m,1),'b');
%     scatter3(reshape(TE_Right_proj(:,:,1),n*m,1),reshape(TE_Right_proj(:,:,2),n*m,1),reshape(TE_Right_proj(:,:,3),n*m,1),'mx');
%     scatter3(reshape(TE_Left_proj(:,:,1),n*m,1),reshape(TE_Left_proj(:,:,2),n*m,1),reshape(TE_Left_proj(:,:,3),n*m,1),'kx');

    
    %
    %     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),HS_vec(:,:,1),HS_vec(:,:,2),HS_vec(:,:,3));
    %     quiver3(LEmidpoint(:,:,1),LEmidpoint(:,:,2),LEmidpoint(:,:,3),LE_vec(:,:,1),LE_vec(:,:,2),LE_vec(:,:,3));
    % DVE Normal Vector
        quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),DVE_norm(:,:,1),DVE_norm(:,:,2),DVE_norm(:,:,3),0,'g');
    
%     quiver3(LE_Left(:,:,1),LE_Left(:,:,2),LE_Left(:,:,3),LE_vec(:,:,1),LE_vec(:,:,2),LE_vec(:,:,3),0);
%     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),xsi_vec(:,:,1),xsi_vec(:,:,2),xsi_vec(:,:,3),0);
%     quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),DVE_norm(:,:,1),DVE_norm(:,:,2),DVE_norm(:,:,3));
%     quiver3(LE_Left(:,:,1),LE_Left(:,:,2),LE_Left(:,:,3),LE_vec_local(:,:,1),LE_vec_local(:,:,2),LE_vec_local(:,:,3),0);

    %     quiver3(TE_Left_proj(:,:,1),TE_Left_proj(:,:,2),TE_Left_proj(:,:,3),TE_vec_proj(:,:,1),TE_vec_proj(:,:,2),TE_vec_proj(:,:,3));
    
%         plot3(panel(:,1),panel(:,2),panel(:,3))
%         surf(panelX,panelY,panelZ)


for i = 1:m
    for j = 1:n
        %computing left-leading edge point in local ref. frame
        tempA(1) = -xsi(i,j) - eta(i,j)*tand(phi_LE(i,j));
        tempA(2) = -eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x1 = tempAA+reshape(CP(i,j,:),1,3,1);
        
        % 		computing left-trailing edge point in local ref. frame
        tempA(1) = xsi(i,j) - eta(i,j)*tand(phi_TE(i,j));
        tempA(2) = -eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x2 = tempAA+reshape(CP(i,j,:),1,3,1);
        
        %computing right-trailing edge point in local ref. frame
        tempA(1) = xsi(i,j) + eta(i,j)*tand(phi_TE(i,j));
        tempA(2) = eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x3 = tempAA+reshape(CP(i,j,:),1,3,1);
        
        %computing right-leading edge point in local ref. frame
        tempA(1) = -xsi(i,j) + eta(i,j)*tand(phi_LE(i,j));
        tempA(2) = eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x4 = tempAA+reshape(CP(i,j,:),1,3,1);
        
        fillX = [x1(1) x2(1) x3(1) x4(1) x1(1)];
        fillY = [x1(2) x2(2) x3(2) x4(2) x1(2)];
        fillZ = [x1(3) x2(3) x3(3) x4(3) x1(3)];
        
        fill3(fillX,fillY,fillZ,'b','FaceAlpha',0.25,'EdgeColor','b')
    end
end
clear x1 x2 x3 x4 tempA tempAA


    hold off
    axis equal
    grid on
end




%%
if debug~=1
    clear chordX chordY chordZ panelX panelY panelZ debug
    clear LE_norm
    clear X1 X2 Y1 Y2 Z1 Z2 AX1 AZ1 AX2 AZ2 C1 C2 e2 e1 chordbase
    clear panel paneldata i j M N
    clear  halfchord halfspan
    
end
