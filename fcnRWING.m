function [vecR] = fcnRWING(valTIMESTEP, SURF, WAKE, FLAG)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(SURF.valNELE*3,1);

len = length(SURF.matCENTER(:,1));

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(SURF.matUINF, SURF.matDVENORM,2);    
else
    [w_wake] = fcnWDVEVEL(SURF.matCENTER, valTIMESTEP, WAKE, SURF, FLAG);

    % Including the wake-induced velocities,
    vecR(end-(len-1):end) = (4*pi).*dot(SURF.matUINF+w_wake, SURF.matDVENORM,2);  

end

end

