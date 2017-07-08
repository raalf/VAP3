function [vecR] = fcnRWING_OL(valNELE, valTIMESTEP, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
    matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVETESWP, vecSYM, valWSIZE, vecDVELE, matVLST, matDVE)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valNELE*5,1);

% len = length(matCENTER(:,1));
% normals = matDVENORM;
% uinf = matUINF;

len = length(matCENTER(:,1)) + length(nonzeros(vecDVELE > 0));
normals = [matDVENORM; matDVENORM(vecDVELE > 0,:)];
uinf = [matUINF; matUINF(vecDVELE > 0,:)];

fpg_le = (matVLST(matDVE(vecDVELE > 0,1),:) + matVLST(matDVE(vecDVELE > 0,2),:))./2;

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(uinf, normals,2);    
else
    [w_wake] = fcnWDVEVEL(matCENTER, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, ...
        vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, valTIMESTEP);
    
    w_wake_le = fcnWDVEVEL(fpg_le, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, ...
        vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, valTIMESTEP);
    
    w_wake = [w_wake; w_wake_le];
    % Including the wake-induced velocities,
    vecR(end-(len-1):end) = (4*pi).*dot(uinf+w_wake, normals,2);  

end

end

