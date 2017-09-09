function [matESHEET] = fcnESHEETS_OL(valNELE, matVLST, matDVE, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW)
% This function is unique to OL, it finds which sheets must be added or subtracted to form the
% finite vortex sheet which induces in the streamwise direction
% matESHEETS - Level 1 columns are edge numbers (1 - 4) and either 1, 0 or -1 based on if the edge sheet
% is added, ignored (90 degree sweep), or subtracted
% Level 2 is phi

matESHEET = zeros(valNELE,4,2);
DBL_EPS = 1e-14;

%% First edge
temp = zeros(valNELE,3,2);
temp(:,:,1) = fcnGLOBSTAR(matVLST(matDVE(:,1),:) - matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
temp(:,:,2) = fcnGLOBSTAR(matVLST(matDVE(:,2),:) - matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
ab = (temp(:,:,2) - temp(:,:,1))./sqrt(sum((temp(:,:,2) - temp(:,:,1)).^2,2));
xdir = repmat([1 0 0], size(ab,1),1);

idx = temp(:,1,1) - temp(:,1,2) >= DBL_EPS;
matESHEET(idx,1,1) = -1;
matESHEET(idx,1,2) = abs(acos(dot(ab(idx), xdir(idx), 2)));

idx = temp(:,1,1) - temp(:,1,2) <= -DBL_EPS;
matESHEET(idx,1,1) = 1;
matESHEET(idx,1,2) = -abs(acos(dot(ab(idx), xdir(idx), 2)));

%% Second edge
idx = sqrt(sum((matVLST(matDVE(:,2),:) - matVLST(matDVE(:,3),:)).^2,2)) > DBL_EPS;
matESHEET(idx,2,1) = 1;

%% Third edge
temp = zeros(valNELE,3,2);
temp(:,:,1) = fcnGLOBSTAR(matVLST(matDVE(:,3),:) - matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
temp(:,:,2) = fcnGLOBSTAR(matVLST(matDVE(:,4),:) - matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
ab = (temp(:,:,2) - temp(:,:,1))./sqrt(sum((temp(:,:,2) - temp(:,:,1)).^2,2));

idx = temp(:,1,1) - temp(:,1,2) >= DBL_EPS;
matESHEET(idx,3,1) = -1;
matESHEET(idx,3,2) = -abs(acos(dot(ab(idx), xdir(idx), 2)));

idx = temp(:,1,1) - temp(:,1,2) <= -DBL_EPS;
matESHEET(idx,3,1) = 1;
matESHEET(idx,3,2) = abs(acos(dot(ab(idx), xdir(idx), 2)));

%% Fourth edge
idx = sqrt(sum((matVLST(matDVE(:,1),:) - matVLST(matDVE(:,4),:)).^2,2)) > DBL_EPS;
matESHEET(idx,4,1) = -1;

end

