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

a1 = zeros(len*2,3);
b1 = zeros(len*2,3);
c1 = zeros(len*2,3);

a2 = zeros(len*2,3);
b2 = zeros(len*2,3);
c2 = zeros(len*2,3);

endpoints = zeros(len,3,4);

%% Leading Edge

% Leading edge coordinates, used to find the midpoint of the leading edge
% A vector is then made between this midpoint and the field point, which is
% then rotated into the DVE reference frame.
endpoints(:,:,1) = matVLST(matDVE(dvenum,1),:); % Left leading edge point
endpoints(:,:,2) = matVLST(matDVE(dvenum,2),:); % Right leading edge point

endpoints(:,:,3) = matVLST(matDVE(dvenum,4),:); % Left trailing edge point
endpoints(:,:,4) = matVLST(matDVE(dvenum,3),:); % Right trailing edge point

xsiA = fcnGLOBSTAR([fpg - mean(endpoints(:,:,1:2),3); fpg - mean(endpoints(:,:,3:4),3)], repmat(vecDVEROLL(dvenum),2,1), repmat(vecDVEPITCH(dvenum),2,1), repmat(vecDVEYAW(dvenum),2,1));

% Bound vortex on the leading edge
idx1 = dvetype == 0 | dvetype == 2 | dvetype == -3 | dvetype == -4;
% Bound vortex at the trailing edge
idx2 = dvetype == 0 | dvetype == -2;
idx12 = [idx1; idx2];

[a1(idx12,:), b1(idx12,:), c1(idx12,:)] = fcnBOUNDIND([vecDVEHVSPN(dvenum(idx1)); vecDVEHVSPN(dvenum(idx2))], [vecDVELESWP(dvenum(idx1)); vecDVETESWP(dvenum(idx2))], xsiA(idx12,:));

% Vortex sheet at the leading and trailing edge
idx3 = true(len,1);
idx4 = dvetype ~= 3 | dvetype ~= -3;
idx34 = [idx3; idx4];
[a2(idx34,:), b2(idx34,:), c2(idx34,:)] = fcnVSIND([vecDVEHVSPN(dvenum(idx3)); vecDVEHVSPN(dvenum(idx4))], [vecDVELESWP(dvenum(idx3)); vecDVETESWP(dvenum(idx4))], xsiA(idx34,:), [vecK(dvenum(idx3)); vecK(dvenum(idx4))]);

%% Summing together the influences from the sheets and filaments

% a3xi = a1(1:len,:) + a2(1:len,:) - a1(len+1:end,:) - a2(len+1:end,:);
% b3xi = b1(1:len,:) + b2(1:len,:) - b1(len+1:end,:) - b2(len+1:end,:);
% c3xi = c1(1:len,:) + c2(1:len,:) - c1(len+1:end,:) - c2(len+1:end,:);

temp = fcnSTARGLOB([a1(1:len,:) + a2(1:len,:) - a1(len+1:end,:) - a2(len+1:end,:);...
    b1(1:len,:) + b2(1:len,:) - b1(len+1:end,:) - b2(len+1:end,:); ...
    c1(1:len,:) + c2(1:len,:) - c1(len+1:end,:) - c2(len+1:end,:)], repmat(vecDVEROLL(dvenum),3,1), repmat(vecDVEPITCH(dvenum),3,1), repmat(vecDVEYAW(dvenum),3,1));

a = temp(1:len,:);
b = temp(len+1:2*len,:);
c = temp((2*len)+1:end,:);

% a = fcnSTARGLOB(a3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
% b = fcnSTARGLOB(b3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
% c = fcnSTARGLOB(c3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

end



















