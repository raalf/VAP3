function [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, vecUINF)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valNELE*3,1);

if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    len = length(matCENTER(:,1));
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1), matDVENORM,2);    
else
    % some shit
end

end

