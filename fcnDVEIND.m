function [a, b, c] = fcnDVEIND(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP)

% This function gives the influence of a DVE.

% DVE_type = 0 DVE has vortex filaments at leading and
%			 	trailing edge, usually lifting surface DVE
% DVE_type = 1 DVE has no vortex filaments at leading and
% 				trailing edge, usually wake DVE
% DVE_type = 2 DVE has vortex filament at leading edge,
% 				but not at trailing edge
% DVE_type =-2 DVE has vortex filament at trailing edge,
% 				but not at leading edge
% DVE_type = 3 DVE is a semi infinite vortex sheet without a
% 				vortex filaments at leading and trailing edge
% DVE_type =-3 DVE is a semi infinite vortex sheet with a
% 				vortex filaments at its leading edge

% INPUT:
%   dvenum - vector of influencing DVEs
%   dvetype - type of DVE, corresponds to above dvenum
%   fpg - list of field points we influence on, length is same as dvenum
%   vecK - singularity factor of the DVEs
%   The rest can be found in the Variables.txt if you are unfamiliar
% OUTPUT:
%   a,b,c - influence coefficients (each are (x,y,z))

% Vector to field point from DVE control point
% rA = fpg - matCENTER(dvenum,:);

% Rotate the above vector into local reference frame
% xsiA = fcnGLOBSTAR(rA, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

len = length(dvenum);

a1le = zeros(len,3);
b1le = zeros(len,3);
c1le = zeros(len,3);

a1te = zeros(len,3);
b1te = zeros(len,3);
c1te = zeros(len,3);

endpoints = zeros(len,3,2);

%% Leading Edge

% Leading edge coordinates, used to find the midpoint of the leading edge
% A vector is then made between this midpoint and the field point, which is
% then rotated into the DVE reference frame.
endpoints(:,:,1) = matVLST(matDVE(dvenum,1),:); % Left leading edge point
endpoints(:,:,2) = matVLST(matDVE(dvenum,2),:); % Right leading edge point
xsiA = fcnGLOBSTAR(fpg - mean(endpoints,3), vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

% Bound vortex on the leading edge
idx1 = dvetype == 0 | dvetype == 2 | dvetype == -3 | dvetype == -4;
[a1le(idx1,:), b1le(idx1,:), c1le(idx1,:)] = fcnBOUNDIND(vecDVEHVSPN(dvenum(idx1)), vecDVELESWP(dvenum(idx1)), xsiA(idx1,:));

% Vortex sheet at leading edge
[a2le, b2le, c2le] = fcnVSIND(vecDVEHVSPN(dvenum), vecDVELESWP(dvenum), xsiA, vecK(dvenum));

clear endpoints

%% Trailing Edge

% Trailing edge coordinates, used to find the midpoint of the trailing edge
% A vector is then made between this midpoint and the field point, which is
% then rotated into the DVE reference frame.
endpoints(:,:,1) = matVLST(matDVE(dvenum,4),:); % Left leading edge point
endpoints(:,:,2) = matVLST(matDVE(dvenum,3),:); % Right leading edge point
xsiA = fcnGLOBSTAR(fpg - mean(endpoints,3), vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

% Bound vortex at the trailing edge
idx2 = dvetype == 0 | dvetype == -2;
[a1te(idx2,:), b1te(idx2,:), c1te(idx2,:)] = fcnBOUNDIND(vecDVEHVSPN(dvenum(idx2)), vecDVETESWP(dvenum(idx2)), xsiA(idx2,:));

% Vortex sheet at the trailing edge
idx3 = dvetype ~= 3 | dvetype ~= -3;
[a2te, b2te, c2te] = fcnVSIND(vecDVEHVSPN(dvenum(idx3)), vecDVETESWP(dvenum(idx3)), xsiA(idx3,:), vecK(dvenum(idx3)));

%% Summing together the influences from the sheets and filaments

a3xi = a1le + a2le - a1te - a2te;
b3xi = b1le + b2le - b1te - b2te;
c3xi = c1le + c2le - c1te - c2te;

a = fcnSTARGLOB(a3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
b = fcnSTARGLOB(b3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
c = fcnSTARGLOB(c3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

end



















