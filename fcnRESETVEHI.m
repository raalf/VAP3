function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)
% Function to reset the vehicle location back to an origin location of
% [0,0,0]

trans = INPU.matVEHORIG - [0,0,0]; % Vector of how much to translate all points


SURF.matEALST = SURF.matEALST - trans;
SURF.matCGLST = SURF.matCGLST - trans;
SURF.vecWINGCG = SURF.vecWINGCG - trans;
VEHI.vecFUSECG = VEHI.vecFUSECG - trans;
VEHI.vecWINGCG = VEHI.vecWINGCG - trans;
VEHI.vecPAYLCG = VEHI.vecPAYLCG - trans;
VEHI.vecPROPLOC = VEHI.vecPROPLOC - trans;
SURF.matAEROCNTR = SURF.matAEROCNTR - trans;

SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) = SURF.matTRIMORIG(FLAG.vecTRIMABLE == 1,:) - trans;

SURF.matCENTER = SURF.matCENTER - trans;
SURF.matVLST = SURF.matVLST - trans;
SURF.matNPVLST = SURF.matNPVLST - trans;
INPU.vecVEHCG = INPU.vecVEHCG - trans;
INPU.matVEHORIG = INPU.matVEHORIG - trans;

end