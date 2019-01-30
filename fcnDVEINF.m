function [a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVESYM)
% This function gives the influence of a DVE, accounting for symmetry.

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
%   a,b,c - influence coefficients (each are (x,y,z)), accounting for symmetry

% T.D.K 2016-09-28 ROTHWELL STREET, AURORA, ONTARIO, CANADA L4G-0V8

[a, b, c] = fcnDVEIND(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP);

if any(vecDVESYM) == 1    
    idx = vecDVESYM(dvenum);
    [as, bs, cs] = fcnDVEIND(dvenum(idx), dvetype(idx), fpg(idx,:), vecK, matDVE, [matVLST(:,1) -matVLST(:,2) matVLST(:,3)], ...
                            vecDVEHVSPN, vecDVEHVCRD,-vecDVEROLL, vecDVEPITCH, -vecDVEYAW, -vecDVELESWP, -vecDVETESWP);
    
    a(idx,:) = a(idx,:) + as;
    b(idx,:) = b(idx,:) - bs;
    c(idx,:) = c(idx,:) + cs; 
end

end

