function [OUTP] = fcnSTABILITYFORCES(INPU, COND, OUTP)
% Function to find the forces and moments on the stability axes and turn
% them into their dimensional and non-dimensional derivative form.

qS = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF*INPU.vecAREA;

% Reference directions for stability axes. X is negative since VAP moves
% the vehicle in the negative X direction
xref = [-1 0 0];
zref = [0 0 -1];

liftref = [-sind(COND.vecVEHALPHA) 0 cosd(COND.vecVEHALPHA)];
dragref = [cosd(COND.vecVEHALPHA) 0 sind(COND.vecVEHALPHA)];

% X, Z forces acting on vehicle
OUTP.vecVEHXFORCE = dot(qS*OUTP.vecCL(end).*liftref,xref) + dot(qS*OUTP.vecCDI(end).*dragref,xref);
OUTP.vecVEHZFORCE = dot(qS*OUTP.vecCL(end).*liftref,zref) + dot(qS*OUTP.vecCDI(end).*dragref,zref);

