function [a, b, c, d, e] = fcnDVEIND_OL(dvenum_all, dvetype_all, fpg_all, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, matESHEETS)

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

chunk_sz = 4e6;

num_pts = length(dvenum_all);

a = zeros(num_pts,3);
b = a;
c = a;
d = a;
e = a;

for i = 1:chunk_sz:num_pts
    
    idx_chunk = [i:i + chunk_sz - 1];
    
    if idx_chunk(end) > num_pts
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
    
    b2te = zeros(len,3);
    c2te = zeros(len,3);
    d2r = zeros(len,3);
    e2r = zeros(len,3);
    d2l = zeros(len,3);
    e2l = zeros(len,3);
    
    %% Leading Edge    
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,1),:) + matVLST(matDVE(dvenum,2),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    [~, b2le, c2le] = fcnVSIND(vecDVEHVSPN(dvenum), vecDVEHVCRD(dvenum), vecDVELESWP(dvenum), xsiA, vecK(dvenum));
    
    %% Trailing Edge
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,3),:) + matVLST(matDVE(dvenum,4),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    idx3 = dvetype ~= 3 & dvetype ~= -3;
    if any(idx3)
        [~, b2te(idx3,:), c2te(idx3,:)] = fcnVSIND(vecDVEHVSPN(dvenum(idx3)), vecDVEHVCRD(dvenum(idx3)), vecDVETESWP(dvenum(idx3)), xsiA(idx3,:), vecK(dvenum(idx3)));
    end
    
    %% Right to left vortex sheets for HDVEs
    idx4 = dvetype == 0;
    
    % Edge 1
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,1),:) + matVLST(matDVE(dvenum,2),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    hspan = fcnGLOBSTAR((matVLST(matDVE(dvenum(idx4),1),:) - matVLST(matDVE(dvenum(idx4),2),:)), vecDVEROLL(dvenum(idx4)), vecDVEPITCH(dvenum(idx4)), vecDVEYAW(dvenum(idx4)));
    [~, d2f, e2f] = fcnVSIND(abs(hspan(:,1))./2, vecDVEHVSPN(dvenum(idx4)), matESHEETS(dvenum(idx4),1, 2), [-xsiA(idx4,2) xsiA(idx4,1) xsiA(idx4,3)], vecK(dvenum(idx4)));
    d2f = [d2f(:,2) -d2f(:,1) d2f(:,3)];
    e2f = [e2f(:,2) -e2f(:,1) e2f(:,3)];
    
    % Edge 2
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,2),:) + matVLST(matDVE(dvenum,3),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    hspan = fcnGLOBSTAR((matVLST(matDVE(dvenum(idx4),3),:) - matVLST(matDVE(dvenum(idx4),2),:)), vecDVEROLL(dvenum(idx4)), vecDVEPITCH(dvenum(idx4)), vecDVEYAW(dvenum(idx4)));
    [~, d2r, e2r] = fcnVSIND(abs(hspan(:,1))./2, vecDVEHVSPN(dvenum(idx4)), matESHEETS(dvenum(idx4),2, 2), [-xsiA(idx4,2) xsiA(idx4,1) xsiA(idx4,3)], vecK(dvenum(idx4)));
    d2r = [d2r(:,2) -d2r(:,1) d2r(:,3)];
    e2r = [e2r(:,2) -e2r(:,1) e2r(:,3)];
    
    % Edge 3
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,3),:) + matVLST(matDVE(dvenum,4),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    hspan = fcnGLOBSTAR((matVLST(matDVE(dvenum(idx4),3),:) - matVLST(matDVE(dvenum(idx4),4),:)), vecDVEROLL(dvenum(idx4)), vecDVEPITCH(dvenum(idx4)), vecDVEYAW(dvenum(idx4)));   
    [~, d2re, e2re] = fcnVSIND(abs(hspan(:,1))./2, vecDVEHVSPN(dvenum(idx4)), matESHEETS(dvenum(idx4),3, 2), [-xsiA(idx4,2) xsiA(idx4,1) xsiA(idx4,3)], vecK(dvenum(idx4)));
    d2re = [d2re(:,2) -d2re(:,1) d2re(:,3)];
    e2re = [e2re(:,2) -e2re(:,1) e2re(:,3)]; 
    
    % Edge 4
    xsiA = fcnGLOBSTAR(fpg - (matVLST(matDVE(dvenum,1),:) + matVLST(matDVE(dvenum,4),:))./2, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    hspan = fcnGLOBSTAR((matVLST(matDVE(dvenum(idx4),4),:) - matVLST(matDVE(dvenum(idx4),1),:)), vecDVEROLL(dvenum(idx4)), vecDVEPITCH(dvenum(idx4)), vecDVEYAW(dvenum(idx4)));
    [~, d2l, e2l] = fcnVSIND(abs(hspan(:,1))./2, vecDVEHVSPN(dvenum(idx4)), matESHEETS(dvenum(idx4),4, 2), [-xsiA(idx4,2) xsiA(idx4,1) xsiA(idx4,3)], vecK(dvenum(idx4)));
    d2l = [d2l(:,2) -d2l(:,1) d2l(:,3)];
    e2l = [e2l(:,2) -e2l(:,1) e2l(:,3)];

    if isempty(e2r); e2r = zeros(size(c1le)); end
    if isempty(e2l); e2l = zeros(size(c1le)); end
    if isempty(d2r); d2r = zeros(size(c1le)); end
    if isempty(d2l); d2l = zeros(size(c1le)); end   
    if isempty(e2re); e2re = zeros(size(c1le)); end
    if isempty(e2f); e2f = zeros(size(c1le)); end
    if isempty(d2re); d2re = zeros(size(c1le)); end
    if isempty(d2f); d2f = zeros(size(c1le)); end 
    
    %% Summing together the influences from the sheets and filaments
    a3xi = a1le - a1te;
    b3xi = b1le + b2le - b1te - b2te;
    c3xi = c1le + c2le - c1te - c2te;
    d3xi = matESHEETS(dvenum,1,1).*d2f + matESHEETS(dvenum,2,1).*d2r + matESHEETS(dvenum,3,1).*d2re + matESHEETS(dvenum,4,1).*d2l;
    e3xi = matESHEETS(dvenum,1,1).*e2f + matESHEETS(dvenum,2,1).*e2r + matESHEETS(dvenum,3,1).*e2re + matESHEETS(dvenum,4,1).*e2l;
   
    a(idx_chunk,:) = fcnSTARGLOB(a3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    b(idx_chunk,:) = fcnSTARGLOB(b3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    c(idx_chunk,:) = fcnSTARGLOB(c3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    d(idx_chunk,:) = fcnSTARGLOB(d3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    e(idx_chunk,:) = fcnSTARGLOB(e3xi, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));
    
    
end

end
