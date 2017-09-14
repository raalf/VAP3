function [a, b, c, d, e, f] = fcnDVEINF_OL(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM, matESHEETS)
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

% This function does not account for symmetry well, it is all or nothing with symmetry,
% but it really should be wing-by-wing

% T.D.K 2016-09-28 ROTHWELL STREET, AURORA, ONTARIO, CANADA L4G-0V8

[a, b, c, d, e, f] = fcnDVEIND_OL(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, matESHEETS);

if any(vecSYM) == 1
    
    matVLSTs = [matVLST(:,1) -matVLST(:,2) matVLST(:,3)];
    
    [as, bs, cs, ds, es, fs] = fcnDVEIND_OL(dvenum, dvetype, fpg, vecK, matDVE, matVLSTs, vecDVEHVSPN, vecDVEHVCRD,-vecDVEROLL, vecDVEPITCH, -vecDVEYAW, -vecDVELESWP, -vecDVETESWP, matESHEETS);
    
    a = a + as;
    b = b - bs;
    c = c + cs;
    d = d - ds;
    e = e + es;
    f = f + fs;
    
end

end

