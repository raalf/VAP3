function [OUTP] = fcnBEAMDISPLACE(OUTP, INPU, SURF, COND, VEHI, FLAG, valTIMESTEP, tempTIME)
% This function computes the spanwise deflection and twist using an
% explicit finite difference method given a loading and structural
% distribution.
%
% INPUTS:
%
% valTIMESTEP - Current timestep number
%
% OUTP.vecLIFTDIST - 1 x n vector with the lift values at each node, where n is
% the number of spanwise stations
%
% OUTP.vecMOMDIST - 1 x n vector with aerodynamic moment at each node, where n
% is the number of spanwise stations
%
% vecSPANAREA - 1 x n vector containing the structural cross sectional area
% at each spanwise location, where n is the number of spanwise stations
%
% INPU.matEIx - n x 3 matrix containing the bending stiffness at each spanwise
% station, where n is the number of spanwise stations. The first column
% represents EIx, the second column EIx', and the third column EIx''
%
% INPU.matGJt - n x 2 matrix containing the torsional stiffness distribution at
% each spanwise node, where n is the number of spanwise stations. The first
% row is GJt and the secon GJt'
%
% vecTORSIONRIGIDITY - 1 x n vector containing the torsional rigidity (GJ)
% at each spanwise station, where n is the number of spanwise stations.
%
% vecMASS2SHEAR - 1 x n vector containing the distances (in m) between the
% center of mass and shear center at each spanwise station, where n is the
% number of spanwise stations
%
% valYMODULUS - Young's modulus of the wing structure in Pa
%
% INPU.valNSELE - Number of spanwise elements
%

[ledves, ~, ~] = find(SURF.vecDVELE(SURF.idxFLEX) > 0);

valDY = 0.5*INPU.vecSPAN/(INPU.valNSELE-1);

temp_y = (0:valDY:0.5*INPU.vecSPAN)';

valSTRUCTDELTIME = COND.valSDELTIME;

%% Interpolate values at structural grid points
if tempTIME == 1
    [SURF, OUTP, COND, INPU] = fcnBEAMFORCE(SURF, OUTP, COND, INPU, VEHI, FLAG, valDY, temp_y, valTIMESTEP);
%     OUTP.vecLIFTDIST_STRUCT = OUTP.vecLIFTDIST./valDY;
%     OUTP.vecMOMDIST_STRUCT = OUTP.vecMOMDIST_NEW./valDY;
end

OUTP.vecDEF = zeros(1,INPU.valNSELE+4);
OUTP.vecTWIST = zeros(1,INPU.valNSELE+4);

valSTRUCTTIME = tempTIME + 2;

%% Beam boundary conditions
% if tempTIME == 1
%     OUTP.matDEF(1:valSTRUCTTIME-1,:) = OUTP.matDEF((OUTP.valSTRUCTITER-1):OUTP.valSTRUCTITER,:);
%     OUTP.matTWIST(1:valSTRUCTTIME-1,:) = OUTP.matTWIST((OUTP.valSTRUCTITER-1):OUTP.valSTRUCTITER,:);
% end

OUTP.vecDEF(3) = SURF.matBEAMLOC(1,3,valTIMESTEP-1); % Zero deflection at root BC
OUTP.vecTWIST(3) = COND.vecVEHPITCH; % Zero twist at root BC

% Assemble load matrix
matLOAD = [OUTP.vecBEAMFORCE, OUTP.vecBEAMMOM];

for yy = 4:(INPU.valNSELE+2)

    %% Geometric property assembly

    % Assemble mass matrix
    matMASS = [INPU.vecLM(yy-2), -INPU.vecLM(yy-2).*SURF.vecLSM(yy-2); -INPU.vecLM(yy-2).*SURF.vecLSM(yy-2), INPU.vecJT(yy-2)];

    % Assemble stiffness matrices
    matK_1 = [INPU.matEIx(yy-2,3), 0; 0, 0];
    matK_2 = [INPU.matEIx(yy-2,2), 0; 0, -INPU.matGJt(yy-2,2)];
    matK_3 = [INPU.matEIx(yy-2,1), 0; 0, -INPU.matGJt(yy-2,1)];
    matB = [5 0; 0 5];

    %% Finite difference relations for partial derivatives

    % Finite difference relations for partial derivatives w.r.t
    % time
    valUDOT = (OUTP.matDEF(valSTRUCTTIME-1,yy) - OUTP.matDEF(valSTRUCTTIME - 2, yy))./valSTRUCTDELTIME;
    valTDOT = (OUTP.matTWIST(valSTRUCTTIME-1,yy) - OUTP.matTWIST(valSTRUCTTIME - 2,yy))./valSTRUCTDELTIME;

    % Finite difference relations for partial derivative of deflection w.r.t Y
    valU_yy = (OUTP.matDEF(valSTRUCTTIME-1,yy+1) - 2*OUTP.matDEF(valSTRUCTTIME-1,yy) + OUTP.matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^2;
    valU_yyy = (OUTP.matDEF(valSTRUCTTIME-1,yy+2) - 3*OUTP.matDEF(valSTRUCTTIME-1,yy+1) + 3*OUTP.matDEF(valSTRUCTTIME-1,yy)- ...
        OUTP.matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^3;
    valU_yyyy = (OUTP.matDEF(valSTRUCTTIME-1,yy+2) - 4*OUTP.matDEF(valSTRUCTTIME-1,yy+1) + 6*OUTP.matDEF(valSTRUCTTIME-1,yy) - ...
        4*OUTP.matDEF(valSTRUCTTIME-1,yy-1) + OUTP.matDEF(valSTRUCTTIME-1,yy-2))/(valDY)^4;

    % Finite difference relations for partial derivative of twist w.r.t Y
    valTHETA_y = (OUTP.matTWIST(valSTRUCTTIME-1,yy+1) - OUTP.matTWIST(valSTRUCTTIME-1,yy-1))/(2*valDY);
    valTHETA_yy = (OUTP.matTWIST(valSTRUCTTIME-1,yy+1) - 2*OUTP.matTWIST(valSTRUCTTIME-1,yy) + OUTP.matTWIST(valSTRUCTTIME-1,yy-1))/(valDY^2);

    %% Solve matrix equation

    % Temp variable with the wing deflection and twist stored as a matrix. The
    % first row is the deflection, w/ each column as a spanwise station. The
    % second row is the twist, w/ each column as a spanwise station.

    tempTWISTBEND = 2.*[OUTP.matDEF(valSTRUCTTIME-1,yy); OUTP.matTWIST(valSTRUCTTIME-1,yy)] - [OUTP.matDEF(valSTRUCTTIME-2,yy); OUTP.matTWIST(valSTRUCTTIME-2,yy)] ...
        + (valSTRUCTDELTIME^2).*inv(matMASS)*([matLOAD(yy-2,1); matLOAD(yy-2,2)] - matK_1*[valU_yy; 0] - matK_2*[valU_yyy; valTHETA_y] -...
        matK_3*[valU_yyyy; valTHETA_yy] - matB*[valUDOT; valTDOT]);

    % Output result of deflection and twist to separate vectors
    OUTP.vecDEF(yy) = tempTWISTBEND(1,:);
    OUTP.vecTWIST(yy) = tempTWISTBEND(2,:);

end

OUTP.vecDEF(INPU.valNSELE+3) = 2*OUTP.vecDEF(INPU.valNSELE+2)...
    -OUTP.vecDEF(INPU.valNSELE+1); % BC for deflection one element beyond wing (positive span direction)

OUTP.vecDEF(INPU.valNSELE+4) = 3*OUTP.vecDEF(INPU.valNSELE+2)...
    -2*OUTP.vecDEF(INPU.valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)

OUTP.vecDEF(2) = OUTP.vecDEF(4); % BC for deflection one element beyond root (negative span direction)

OUTP.vecTWIST(INPU.valNSELE+3) = OUTP.vecTWIST(INPU.valNSELE+1); % BC for twist one element beyond wing tip (positive span direction)

OUTP.matDEF(valSTRUCTTIME,:) = OUTP.vecDEF;
OUTP.matTWIST(valSTRUCTTIME,:) = OUTP.vecTWIST;

% Spanwise deflection and twist wrt structural timestep
OUTP.vecDEF = OUTP.matDEF(end,:);
OUTP.vecTWIST = OUTP.matTWIST(end,:);

% Spanwise deflection and twist wrt to global timestep
OUTP.matDEFGLOB(valTIMESTEP,:) = interp1(temp_y,OUTP.matDEF(valSTRUCTTIME,3:end-2),SURF.vecSPANLOC);
OUTP.matTWISTGLOB(valTIMESTEP,:) = interp1(temp_y,OUTP.matTWIST(valSTRUCTTIME,3:end-2),SURF.vecSPANLOC);

for i = 2:length(SURF.vecSPANLOC)
    OUTP.matSLOPE(valTIMESTEP,i) = asin((OUTP.matDEFGLOB(valTIMESTEP,i)-OUTP.matDEFGLOB(valTIMESTEP,i-1))/(SURF.vecSPANLOC(i)-SURF.vecSPANLOC(i-1)));
end


end