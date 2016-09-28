function [a, b, c] = fcnDVEIND(dvenum, dvetype, fpg, vecK, matDVE, matCENTER, matDVEROLL, matDVEPITCH, matDVELESWP, matDVETESWP, matDVEYAW, matDVEHVSPN, vecDVEHVCRD)

% Input DVE number and field point and get back a,b,c

[aloc, bloc, cloc] = fcnBOUNDIND(endpoints, phi, yaw, fpl)
[aloc, bloc, cloc] = fcnVSIND(endpoints, phi, yaw, fpl, k)

end

