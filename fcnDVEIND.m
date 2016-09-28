function [a, b, c] = fcnDVEIND(dvenum, dvetype, fpg, vecK, matDVE, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecDVEHVCRD)

% Input DVE number and field point and get back a,b,c

% Vector to field point from DVE control point
rA = fpg - matCENTER(dvenum,:);

% Rotate the above vector into local reference frame
xsiA = fcnGLOBSTAR(rA, vecDVEROLL(dvenum), vecDVEPITCH(dvenum), vecDVEYAW(dvenum));


idx1 = dvetype == 0 | dvetype == 2 | dvetype == -3 | dvetype == -4;

endpoints(:,:,1) = matVLST(matDVE(dvenum,1),:); % Left leading edge point
endpoints(:,:,2) = matVLST(matDVE(dvenum,2),:); % Right leading edge point

[a1le(idx1), b1le(idx1), c1le(idx1)] = fcnBOUNDIND(endpoints(idx1,:,:), vecDVELESWP(dvenum & idx1), 0, xsiA(idx1,:));


[aloc, bloc, cloc] = fcnBOUNDIND(endpoints, phi, yaw, fpl)
[aloc, bloc, cloc] = fcnVSIND(endpoints, phi, yaw, fpl, k)

end

