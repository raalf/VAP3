function [panel] = fcnPanelCorners(rLE,tLE,rchord,tchord,repsilon,tepsilon)
%=========================================================================%
%This function creates one panel for DesignTools.

%Uses: paneldata structure for rLE coordinates,rchord,repsilon,tLE
%coordinates,tchord and tepsilon
%Uses: root TE value
%Uses: panelcount value

%Returns: the coordinates of the four corners of the panel in the form:
%X1LE, Y1LE, Z1LE, X2LE, Y2LE, Z2LE, X2TE, Y2TE, Z2TE, X1TE, Y1TE, Z1TE
%=========================================================================%

%find span and dihedral angle of panel (nu)
panelspan = sqrt((tLE(3)-rLE(3))^2 + (tLE(2)-rLE(2))^2);
nu = asin((tLE(3)-rLE(3))/panelspan);

%build panel
% LE root (X,Y,Z), from the paneldata from the design.txt file if allign =
% false
panel(1:3) = [rLE(1),rLE(2),rLE(3)];

% LE tip (X,Y,Z), from the paneldata from the design.txt file
panel(4:6)= [tLE(1),tLE(2),tLE(3)];

% TE tip (X,Y,Z), take LE tip and add
panel(7:9) = panel(4:6) + [tchord*cos(tepsilon),tchord*sin(tepsilon)*sin(nu),tchord*-sin(tepsilon)*cos(nu)];

% TE root (X,Y,Z)

% the TE root location is built off the LE root specified
% in the design.txt with the given root chord, epsilon and the calculated
% dihedral (nu)
panel(10:12) = panel(1:3) + [rchord*cos(repsilon),rchord*sin(repsilon)*sin(nu),rchord*-sin(repsilon)*cos(nu)];

end