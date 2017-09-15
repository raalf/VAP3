function [matESHEET] = fcnESHEETS_OL(valNELE, matVLST, matDVE, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW)
% This function is unique to OL, it finds which sheets must be added or subtracted to form the
% finite vortex sheet which induces in the streamwise direction
% matESHEETS - Level 1 columns are edge numbers (1 - 4) and either 1, 0 or -1 based on if the edge sheet
% is added, ignored (90 degree sweep), or subtracted
% Level 2 is phi

matESHEET = zeros(valNELE,4,2);
DBL_EPS = 1e-14;

%% First edge
temp = fcnGLOBSTAR(matVLST(matDVE(:,2),:) - matVLST(matDVE(:,1),:), vecDVEROLL, vecDVEPITCH, vecDVEYAW);
ab = temp./sqrt(sum(temp.^2,2));

idx = ab(:,1) <= -DBL_EPS;
matESHEET(idx,1,1) = -1;

idx = ab(:,1) >= DBL_EPS;
matESHEET(idx,1,1) = 1;

idx = abs(ab(:,1)) >= DBL_EPS;
matESHEET(idx,1,2) = -atan(ab(idx,1)./ab(idx,2));

%% Second edge
idx = sqrt(sum((matVLST(matDVE(:,2),:) - matVLST(matDVE(:,3),:)).^2,2)) > DBL_EPS;
matESHEET(idx,2,1) = 1;

%% Third edge
temp = fcnGLOBSTAR(matVLST(matDVE(:,3),:) - matVLST(matDVE(:,4),:), vecDVEROLL, vecDVEPITCH, vecDVEYAW);
ab = temp./sqrt(sum(temp.^2,2));

idx = ab(:,1) <= -DBL_EPS;
matESHEET(idx,3,1) = 1;

idx = ab(:,1) >= DBL_EPS;
matESHEET(idx,3,1) = -1;

idx = abs(ab(:,1)) >= DBL_EPS;
matESHEET(idx,3,2) = -atan(ab(idx,1)./ab(idx,2));

%% Fourth edge
idx = sqrt(sum((matVLST(matDVE(:,1),:) - matVLST(matDVE(:,4),:)).^2,2)) > DBL_EPS;
matESHEET(idx,4,1) = -1;

end

