%% INPUTS

%V0 - before fixing spanwise interp

clc
clear
clf
X1 = 0;
X2 = 0;
Z1 = 0;
Z2 = 2.3;
Y1 = 0;
Y2 = 0.01;
C1 = 2;
C2 = 1;
e2 = 45;


paneldata.rLE = [X1,Y1,Z1];
paneldata.tLE = [X2,Y2,Z2];
paneldata.rchord = C1;
paneldata.tchord = C2;

paneldata.repsilon = 0;
paneldata.tepsilon = deg2rad(e2);

panel = reshape(paneldesign(paneldata),3,4)';
panelX = [panel([1;4],1),panel([2;3],1)];
panelY = [panel([1;4],2),panel([2;3],2)];
panelZ = [panel([1;4],3),panel([2;3],3)];
%
n = 6; %chordwise
m = 4; %spanwise

% NChordwise and NSpanwise
% Generate extra points to find control point
N = n*2;
M = m*2;


% split Root Chord to chordwise elements
% rtX = linspace(


% Vq = interp2(panelX,panelY,panelZ,1:2,1:2);

% Root Airfoil Points
AX1 = linspace(panel(1,1),panel(4,1),M+1)';
AZ1 = linspace(panel(1,3),panel(4,3),M+1)';

% Tip Airfoil Points
AX2 = linspace(panel(2,1),panel(3,1),M+1)';
AZ2 = linspace(panel(2,3),panel(3,3),M+1)';


S1 = [AX1,zeros(M+1,1)+Y1,AZ1];
S2 = [AX2,zeros(M+1,1)+Y2,AZ2];
Span = [Y1 Y2];
Span2 = [0:1/N:1].*(Y2-Y1)+Y1;


for i = 1:M+1
    PX = [S1(i,1) S2(i,1)];
    PY = [S1(i,2) S2(i,2)];
    PZ = [S1(i,3) S2(i,3)];
    
    PX2(i,:) = interp1(Span,PX,Span2);
    PY2(i,:) = interp1(Span,PY,Span2);
    PZ2(i,:) = interp1(Span,PZ,Span2);
end
% surf(PX2,PY2,PZ2);
% surf(PX2(1:2:end,1:2:end),PY2(1:2:end,1:2:end),PZ2(1:2:end,1:2:end));
% colormap white

axis equal





% Calculate Control Points
CP(:,:,1) = PX2(2:2:end,2:2:end);
CP(:,:,2) = PY2(2:2:end,2:2:end);
CP(:,:,3) = PZ2(2:2:end,2:2:end);
%
hold on
scatter3(reshape(CP(:,:,1),n*m,1),reshape(CP(:,:,2),n*m,1),reshape(CP(:,:,3),n*m,1),'xg');
hold off

%
clc
RightofCP(:,:,1) = PX2(2:2:end,3:2:end);
RightofCP(:,:,2) = PY2(2:2:end,3:2:end);
RightofCP(:,:,3) = PZ2(2:2:end,3:2:end);

LEmidpoint(:,:,1) = PX2(1:2:end-1,2:2:end);
LEmidpoint(:,:,2) = PY2(1:2:end-1,2:2:end);
LEmidpoint(:,:,3) = PZ2(1:2:end-1,2:2:end);





% Parameters Calculation
HS = elementwise3Ddist(RightofCP,CP);
HC = elementwise3Ddist(LEmidpoint,CP);
Area = HS.*HC.*4;



% Leading Edge Vectors
% Find leading edge vectors across entire panel
LE_vec = [PX2(1:2:end-1,2) PY2(1:2:end-1,2) PZ2(1:2:end-1,2)]-[PX2(1:2:end-1,1) PY2(1:2:end-1,1) PZ2(1:2:end-1,1)];
% Normalize LE vector
LE_norm = ((LE_vec(:,1).^2+LE_vec(:,2).^2+LE_vec(:,3).^2).^0.5);
LE_vec = LE_vec./repmat(LE_norm,1,3);
% Repeat results n times for entire panel
LE_vec = repmat(reshape(LE_vec,m,1,3),1,n);


% Find Halfspan Vector. Vector from Control Point to LE midpoint
HS_vec = [LEmidpoint - CP];
% Normalize Halfspan vector
HS_norm = (HS_vec(:,:,1).^2+HS_vec(:,:,2).^2+HS_vec(:,:,3).^2).^0.5;
HS_norm = repmat(HS_norm,1,1,3);
HS_vec = (HS_vec./HS_norm);

% Find Normal vector of the DVE. 
% Leading Edge vector CROSS Halfspan vector
DVE_norm = cross(LE_vec,HS_vec,3);


% Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
roll = -atand(DVE_norm(:,:,2)./DVE_norm(:,:,3))

% Pitch in Degrees
% arcsin ( X component of DVE normal vector ) 
pitch  = asind(DVE_norm(:,:,1));

% Yaw = 0 Degrees
yaw = zeros(m,n);



%
hold on


quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),HS_vec(:,:,1),HS_vec(:,:,2),HS_vec(:,:,3));
quiver3(LEmidpoint(:,:,1),LEmidpoint(:,:,2),LEmidpoint(:,:,3),LE_vec(:,:,1),LE_vec(:,:,2),LE_vec(:,:,3));
quiver3(CP(:,:,1),CP(:,:,2),CP(:,:,3),DVE_norm(:,:,1),DVE_norm(:,:,2),DVE_norm(:,:,3));



% plot3(panel(:,1),panel(:,2),panel(:,3))


surf(panelX,panelY,panelZ)

hold off











% clear PX PY PZ S1 S2 Span Span2 X1 X2 Y1 Y2 Z1 Z2 AX1 AZ1 AX2 AZ2 C1 C2 e2
% clear panel paneldata i j M N
% clear  halfchord halfspan
