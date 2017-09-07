function [FW] = generateWakeDVEs(FW, Output,Temp,ii,j)
% BASED ON ALTON'S generateDVEs FUNCTION!!! REPURPOSED FOR WAKE GEOM
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

% Function is now called one element at a time! 16/01/2016 (Bill)
showplot =0;       %1 = show plot, 0 = hide plot
debug = 0;          %1 = clear all the variables at the end, Only show final matrices

% Panel = Aircraft.Surface(ii).Panel;

% X1 = Panel.X(1); Y1 = Panel.Y(1); Z1 = Panel.Z(1); C1 = Panel.Chord(1); e1 = Panel.Twist(1);
% X2 = Panel.X(2); Y2 = Panel.Y(2); Z2 = Panel.Z(2); C2 = Panel.Chord(2); e2 = Panel.Twist(2);

n = 1;% making one wake element at a time
m = 1; %for wake elements

%forward index to TE row of elements, then add j, j is the spanwise index
jj = (FW.Panels(ii).n * FW.m) - FW.Panels(ii).n  + j;
% ==============================================================================

% Calculate twist at leading edge with dihedral correctly modelled
% paneldata.rLE = [X1,Y1,Z1];
% paneldata.tLE = [X2,Y2,Z2];
% paneldata.rchord = C1;
% paneldata.tchord = C2;
%
% paneldata.repsilon = deg2rad(e1);
% paneldata.tepsilon = deg2rad(e2);

% Read panel corners
% panel = reshape(paneldesign(paneldata),3,4)';
% panelX = [panel([1;4],1),panel([2;3],1)];
% panelY = [panel([1;4],2),panel([2;3],2)];
% panelZ = [panel([1;4],3),panel([2;3],3)];

%element corners. these correspond to the previous timestep panel TEs and the
%new timestep panel TEs. 
% if Temp.timestep == 0 %need to change '-' to 'n' in the name of previous timestep struct
%     panelX = [FW.Panels(ii).TECoord(1,1),FW.Panels(ii).TECoord(2,1);Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(1,1),Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(2,1)];
%     panelY = [FW.Panels(ii).TECoord(1,2),FW.Panels(ii).TECoord(2,2);Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(1,2),Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(2,2)];
%     panelZ = [FW.Panels(ii).TECoord(1,3),FW.Panels(ii).TECoord(2,3);Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(1,3),Output.(Temp.AI).(strrep(strcat('timestep',num2str(Temp.timestep-1)),'-','N'))(ii).TECoord(2,3)];
%     
% else %i
%     panelX = [FW.Panels(ii).TECoord(1,1),FW.Panels(ii).TECoord(2,1);Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(1,1),Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(2,1)];
%     panelY = [FW.Panels(ii).TECoord(1,2),FW.Panels(ii).TECoord(2,2);Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(1,2),Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(2,2)];
%     panelZ = [FW.Panels(ii).TECoord(1,3),FW.Panels(ii).TECoord(2,3);Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(1,3),Output.(Temp.AI).(strcat('timestep',num2str(Temp.timestep-1)))(ii).TECoord(2,3)];
% end


    panelX = [FW.Panels(ii).DVE.TECoordL(jj,1),FW.Panels(ii).DVE.TECoordR(jj,1);Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordL(jj,1),Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordR(jj,1)];
    panelY = [FW.Panels(ii).DVE.TECoordL(jj,2),FW.Panels(ii).DVE.TECoordR(jj,2);Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordL(jj,2),Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordR(jj,2)];
    panelZ = [FW.Panels(ii).DVE.TECoordL(jj,3),FW.Panels(ii).DVE.TECoordR(jj,3);Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordL(jj,3),Output.(Temp.AI)(Temp.timestep+1).TimestepData(ii).DVE.TECoordR(jj,3)];
   
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



%% WAKE DVE Parameters Calculation
% Calculate Control Points, stored in 3D matrix
CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],m,n,3);
CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],m,n,3);
LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],m,n,3);
TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],m,n,3);
TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],m,n,3);

% Filter Leading Edge Point
LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],m,n,3);
LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],m,n,3);

% Left and right side mid-points (used in Relaxed wake)

% xleft = (TE_Left + LE_Left)./2;
% xright = (TE_Right + LE_Right)./2;

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

xleft = (TE_Left_proj + LE_Left)./2;
xright = (TE_Right_proj + LE_Right)./2;


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

%This section gives each wake element a unique identifier, beginning with 1
% and going up along the span, and going up with each timestep. 
name = j;

if ii > 1
    name = name + sum([FW.Panels(1:ii-1).n]);
end

if Temp.timestep > -1
    name = name + (sum([FW.Panels.n])*(Temp.timestep));
end
%ii is panel
%jj is element on panel
indx = Temp.timestep + 1;
FW.Panels(ii).WakeDVE(indx).Index(j) = name;
% FW.Panels(ii).WakeDVE(indx).area = reshape(Area',[],1);
% FW.Panels(ii).WakeDVE(indx).eta = reshape(eta,[],1);
% FW.Panels(ii).WakeDVE(indx).xsi = reshape(xsi,[],1);
% FW.Panels(ii).WakeDVE(indx).roll = reshape(nu,[],1);
% FW.Panels(ii).WakeDVE(indx).pitch = reshape(epsilon,[],1);
% FW.Panels(ii).WakeDVE(indx).yaw = reshape(psi,[],1);
% FW.Panels(ii).WakeDVE(indx).xo = [reshape(CP(:,:,1)',[],1) reshape(CP(:,:,2)',[],1) reshape(CP(:,:,3)',[],1)];
% FW.Panels(ii).WakeDVE(indx).phiLE = reshape(phi_LE',[],1);
% FW.Panels(ii).WakeDVE(indx).phiMID = reshape(phi_MID',[],1);
% FW.Panels(ii).WakeDVE(indx).phiTE = reshape(phi_TE',[],1);
% FW.Panels(ii).WakeDVE(indx).norm = [reshape(DVE_norm(:,:,1)',[],1) reshape(DVE_norm(:,:,2)',[],1) reshape(DVE_norm(:,:,3)',[],1)];

FW.Panels(ii).WakeDVE(indx).singfct(j) = 0;

%% Version2 -> Reformat 2D vector to 1D
count = length(Area(:));
FW.Panels(ii).WakeDVE(indx).area(j) = Area(:);
FW.Panels(ii).WakeDVE(indx).eta(j) = eta(:);
FW.Panels(ii).WakeDVE(indx).xsi(j) = xsi(:);
FW.Panels(ii).WakeDVE(indx).roll(j) = nu(:);
FW.Panels(ii).WakeDVE(indx).pitch(j) = epsilon(:);
FW.Panels(ii).WakeDVE(indx).yaw(j) = psi(:);
FW.Panels(ii).WakeDVE(indx).xo(j,:) = reshape(CP(:),count,3);
FW.Panels(ii).WakeDVE(indx).phiLE(j) = phi_LE(:);
FW.Panels(ii).WakeDVE(indx).phiMID(j) = phi_MID(:);
FW.Panels(ii).WakeDVE(indx).phiTE(j) = phi_TE(:);
FW.Panels(ii).WakeDVE(indx).norm(j,:) = reshape(DVE_norm(:),count,3);
FW.Panels(ii).WakeDVE(indx).xleft(j,:) = reshape(xleft(:), count, 3);
FW.Panels(ii).WakeDVE(indx).xright(j,:) = reshape(xright(:), count, 3);




%% add a b and c to wake DVEs

FW.Panels(ii).WakeDVE(indx).A(j) = FW.Panels(ii).DVE.A(jj);
FW.Panels(ii).WakeDVE(indx).B(j) = FW.Panels(ii).DVE.B(jj);
FW.Panels(ii).WakeDVE(indx).C(j) = FW.Panels(ii).DVE.C(jj);



% Changed by Travis, 2016-04-28
% 503-7 Carlton Street,
% Toronto, Ontario
% Canada, M5B-2M3

% It seems we were using wake elements to calculate the K, instead of the surface elements? For some reason our A, B and C's are still slightly off. I'm looking into that now.

% FW.Panels(ii).WakeDVE(indx).K(j) = FW.Panels(ii).WakeDVE(indx).A(j) + (FW.Panels(ii).WakeDVE(indx).eta(j)*(FW.Panels(ii).WakeDVE(indx).eta(j)/3)*FW.Panels(ii).WakeDVE(indx).C(j));
FW.Panels(ii).WakeDVE(indx).K(j) = FW.Panels(ii).DVE.A(jj) + (FW.Panels(ii).DVE.eta(jj)*(FW.Panels(ii).DVE.eta(jj)/3)*FW.Panels(ii).DVE.C(jj));

% fprintf('K: %f A: %f B: %f C: %f eta: %f\n',FW.Panels(ii).WakeDVE(indx).K(j), FW.Panels(ii).DVE.A(jj), FW.Panels(ii).DVE.B(jj), FW.Panels(ii).DVE.C(jj), FW.Panels(ii).DVE.eta(jj));

%% singularity factor for the wake elements
%the singfct is set equal for all wake elements of a specific wing
%singfct is 1% of the tip element half-span(eta) 

%ii is the panel number we are on

%now we need the tip panel for this wing
for k = 1:size(FW.Joint,2)
    if size(FW.Joint(k).Panel) == 1
        tippanel = FW.Joint(k).Panel;
        %what wing is this the tip panel of?
        wing = FW.Panels(tippanel).Wing;
        %does this match the current panel's wing?
        if FW.Panels(ii).Wing == FW.Panels(tippanel).Wing
        break                  
        end
        
    end
end

tempS = 0.01*FW.Panels(tippanel).DVE.eta(end);
FW.Panels(ii).WakeDVE(indx).singfct(j) = tempS; 
clear i indx tempS


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















