function [OUTP] = fcnELASTICWING(OUTP, INPU, SURF, COND, valTIMESTEP)

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

valDY = 0.5*INPU.vecSPAN/INPU.valNSELE;

temp_y = (0:valDY:0.5*INPU.vecSPAN)';

SURF.vecSPANDIST(end) = temp_y(end);

% INPU.valNSELE = sum(vecN,1)+1;

valSTRUCTDELTIME = COND.valSDELTIME;

% INPU.vecJT = 0.00037078.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST...
%     - 0.01102270.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST...
%     + 0.12838255.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST - 0.73708913.*SURF.vecSPANDIST.*SURF.vecSPANDIST.*SURF.vecSPANDIST...
%     + 2.15067037.*SURF.vecSPANDIST.*SURF.vecSPANDIST - 2.99312818.*SURF.vecSPANDIST + 1.84576176;

%% Interpolate values at structural grid points

[matEIx_interp(:,1)] = linterp(SURF.vecSPANDIST,INPU.matEIx(:,1)',temp_y);
[matEIx_interp(:,2)] = linterp(SURF.vecSPANDIST,INPU.matEIx(:,2)',temp_y);
[matEIx_interp(:,3)] = linterp(SURF.vecSPANDIST,INPU.matEIx(:,3)',temp_y);
[matGJt_interp(:,1)] = linterp(SURF.vecSPANDIST,INPU.matGJt(:,1)',temp_y);
[matGJt_interp(:,2)] = linterp(SURF.vecSPANDIST,INPU.matGJt(:,2)',temp_y);
[INPU.vecLM] = linterp(SURF.vecSPANDIST,INPU.vecLM',temp_y);
[INPU.vecJT] = linterp(SURF.vecSPANDIST,INPU.vecJT',temp_y);
[SURF.vecLSM] = linterp(SURF.vecSPANDIST,SURF.vecLSM',temp_y);
[OUTP.vecLIFTDIST] = linterp(SURF.vecSPANDIST,OUTP.vecLIFTDIST,temp_y);
[OUTP.vecMOMDIST] = linterp(SURF.vecSPANDIST,OUTP.vecMOMDIST',temp_y);

INPU.matEIx = matEIx_interp;
INPU.matGJt = matGJt_interp;

OUTP.vecDEF = zeros(1,INPU.valNSELE+4);
OUTP.vecTWIST = zeros(1,INPU.valNSELE+4);
vecSLOPE = zeros(1,INPU.valNSELE-1);

valSTRUCTTIME = valTIMESTEP;
% valSTRUCTTIME = tempTIME + 2;

%% Beam boundary conditions

% Grab solution from last two time steps 
OUTP.matDEF(valSTRUCTTIME-2:valTIMESTEP-1,:) = OUTP.matDEF_old(end-1:end,:);
OUTP.matTWIST(valSTRUCTTIME-2:valTIMESTEP-1,:) = OUTP.matTWIST_old(end-1:end,:);

% for i=1:(valTIMESTEP-1)
%     OUTP.matDEF(i,3:end-1) = linterp(SURF.vecSPANDIST,OUTP.matDEFGLOB(i,:),temp_y);
%     OUTP.matTWIST(i,3:end-1) = linterp(SURF.vecSPANDIST,OUTP.matTWISTGLOB(i,:),temp_y);
% end

% if tempTIME == 1
%     OUTP.matDEF(1:valSTRUCTTIME-1,:) = OUTP.matDEF((end-1):end,:);
%     OUTP.matTWIST(1:valSTRUCTTIME-1,:) = OUTP.matTWIST((end-1):end,:);
% end

OUTP.vecDEF(3) = 0; % Zero deflection at root BC
OUTP.vecTWIST(3) = 0; % Zero twist at root BC

% Assemble load matrix
matLOAD = [OUTP.vecLIFTDIST' - INPU.vecLM'.*9.81, OUTP.vecMOMDIST' + SURF.vecLSM'.*INPU.vecLM'.*9.81];
% matLOAD = [OUTP.vecLIFTDIST', OUTP.vecMOMDIST'];
% matLOAD(end,:) = [0,0]; 

for yy = 4:(INPU.valNSELE+2)

    %% Geometric property assembly

    % Assemble mass matrix
    matMASS = [INPU.vecLM(yy-2), -INPU.vecLM(yy-2).*SURF.vecLSM(yy-2); -INPU.vecLM(yy-2).*SURF.vecLSM(yy-2), INPU.vecJT(yy-2)];

    % Assemble stiffness matrices
    matK_1 = [INPU.matEIx(yy-2,3), 0; 0, 0];
    matK_2 = [INPU.matEIx(yy-2,2), 0; 0, -INPU.matGJt(yy-2,2)];
    matK_3 = [INPU.matEIx(yy-2,1), 0; 0, -INPU.matGJt(yy-2,1)];
    matB = [0 0; 0 100];

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

    % Calculate angle between DVE and horizontal based on
    % deflection
    vecSLOPE(yy-3) = asin((OUTP.vecDEF(yy)-OUTP.vecDEF(yy-1))/(temp_y(yy-2)-temp_y(yy-3)));

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
OUTP.vecDEF = OUTP.matDEF(valSTRUCTTIME,:);
OUTP.vecTWIST = OUTP.matTWIST(valSTRUCTTIME,:);

vecSLOPE = [0; vecSLOPE'];

% Spanwise deflection and twist wrt to global timestep
OUTP.matDEFGLOB(valTIMESTEP,:) = linterp(temp_y,OUTP.matDEF(valSTRUCTTIME,3:end-1),SURF.vecSPANDIST);
OUTP.matTWISTGLOB(valTIMESTEP,:) = linterp(temp_y,OUTP.matTWIST(valSTRUCTTIME,3:end-1),SURF.vecSPANDIST);

OUTP.matSLOPE(valTIMESTEP,:) = linterp(temp_y(1:end-1),vecSLOPE,SURF.vecSPANDIST(1:end-1));


end