function [a, b, c] = fcnDVEIND_GPU(dvenum_all, dvetype_all, fpg_all, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP)

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

chunk_sz = 18e6;

num_pts = length(dvenum_all);

fp = fopen('memory.txt','a');
fprintf(fp,'%d\n',num_pts);
fclose(fp);

a = zeros(num_pts,3);
b = a;
c = a;

for i = 1:chunk_sz:num_pts
    
    if num_pts <= chunk_sz
        idx_chunk = [1:num_pts];
    elseif i + chunk_sz - 1 <= num_pts
        idx_chunk = [i:i + chunk_sz - 1];
    else
        idx_chunk = [i:num_pts];
    end
    
    fpg = fpg_all(idx_chunk,:);
    dvenum = dvenum_all(idx_chunk,:);
    dvetype = dvetype_all(idx_chunk,:);

    
    
    
    len = length(dvenum);

    a1le = zeros(len,3);
    b1le = zeros(len,3);
    c1le = zeros(len,3);

    a1te = zeros(len,3);
    b1te = zeros(len,3);
    c1te = zeros(len,3);

    b2te = zeros(len,3, 'gpuArray');
    c2te = zeros(len,3, 'gpuArray');

    %% Leading Edge

    % Leading edge coordinates, used to find the midpoint of the leading edge
    % A vector is then made between this midpoint and the field point, which is
    % then rotated into the DVE reference frame.
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,1),:)+matVLST(matDVE(dvenum,2),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

    % xsiA = fcnGLOBSTARGPU(gpuArray(fpg - (matVLST(matDVE(dvenum,1),:)+matVLST(matDVE(dvenum,2),:))./2), gpuArray(vecDVEROLL(dvenum)), gpuArray(vecDVEPITCH(dvenum)), gpuArray(vecDVEYAW(dvenum)));
    % xsiA = gather(xsiA);

    % Bound vortex on the leading edge
    idx1 = dvetype == 0 | dvetype == 2 | dvetype == -3 | dvetype == -4;
    [a1le(idx1,:), b1le(idx1,:), c1le(idx1,:)] = fcnBOUNDIND(vecDVEHVSPN(dvenum(idx1)), vecDVELESWP(dvenum(idx1)), xsiA(idx1,:));

    % Vortex sheet at leading edge
    [~, b2le, c2le] = fcnVSIND(gpuArray(single(vecDVEHVSPN(dvenum))), gpuArray(single(vecDVEHVCRD(dvenum))), gpuArray(single(vecDVELESWP(dvenum))), gpuArray(single(xsiA)), gpuArray(single(vecK(dvenum))), 1); 
    b2le = gather(b2le);
    c2le = gather(c2le);

    %% Trailing Edge

    % Trailing edge coordinates, used to find the midpoint of the trailing edge
    % A vector is then made between this midpoint and the field point, which is
    % then rotated into the DVE reference frame.
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,3),:)+matVLST(matDVE(dvenum,4),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

    % xsiA = fcnGLOBSTARGPU(gpuArray(fpg - (matVLST(matDVE(dvenum,3),:)+matVLST(matDVE(dvenum,4),:))./2), gpuArray(vecDVEROLL(dvenum)), gpuArray(vecDVEPITCH(dvenum)), gpuArray(vecDVEYAW(dvenum)));
    % xsiA = gather(xsiA);

    % Bound vortex at the trailing edge
    idx2 = dvetype == 0 | dvetype == -2;
    [a1te(idx2,:), b1te(idx2,:), c1te(idx2,:)] = fcnBOUNDIND(vecDVEHVSPN(dvenum(idx2)), vecDVETESWP(dvenum(idx2)), xsiA(idx2,:));

    % Vortex sheet at the trailing edge
    idx3 = dvetype ~= 3 & dvetype ~= -3;
    [~, b2te(idx3,:), c2te(idx3,:)] = fcnVSIND(gpuArray(single(vecDVEHVSPN(dvenum(idx3)))), gpuArray(single(vecDVEHVCRD(dvenum(idx3)))), gpuArray(single(vecDVETESWP(dvenum(idx3)))), gpuArray(single(xsiA(idx3,:))), gpuArray(single(vecK(dvenum(idx3)))), 1);
    b2te = gather(b2te);
    c2te = gather(c2te);

    %% Summing together the influences from the sheets and filaments
    a3xi = a1le - a1te;
    b3xi = b1le + b2le - b1te - b2te;
    c3xi = c1le + c2le - c1te - c2te;

    a(idx_chunk,:) = fcnSTARGLOB(a3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    b(idx_chunk,:) = fcnSTARGLOB(b3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    c(idx_chunk,:) = fcnSTARGLOB(c3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));

    % a = fcnSTARGLOBGPU(gpuArray(a3xi), gpuArray(vecDVEROLL(dvenum)), gpuArray(vecDVEPITCH(dvenum)), gpuArray(vecDVEYAW(dvenum)));
    % b = fcnSTARGLOBGPU(gpuArray(b3xi), gpuArray(vecDVEROLL(dvenum)), gpuArray(vecDVEPITCH(dvenum)), gpuArray(vecDVEYAW(dvenum)));
    % c = fcnSTARGLOBGPU(gpuArray(c3xi), gpuArray(vecDVEROLL(dvenum)), gpuArray(vecDVEPITCH(dvenum)), gpuArray(vecDVEYAW(dvenum)));
    % a = gather(a);
    % b = gather(b);
    % c = gather(c);

end