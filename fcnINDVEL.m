function [w_total] = fcnINDVEL(fpg,valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, flagTRI, flagSTEADY, flagGPU)

%returns the induced velocity at a point due to all surface and wake
%elements


% INPUT:
%   pass in a vector of points (num points x 3)
% OUTPUT:
%   w_total - numpoints x 3  induced velocities

w_total = zeros((length(fpg(:,1))),3);
% velocities from surface elements
[w_surf] = fcnSDVEVEL(fpg, valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM, flagGPU);

% velocities from wake elements
[w_wake] = fcnWDVEVEL(fpg, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, valTIMESTEP, flagTRI, flagSTEADY, flagGPU);

% add
w_total = w_surf+w_wake;

