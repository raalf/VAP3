function [matUINF, gust_vel, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valUINF,valGUSTSTART,fpg,gust_vel_old, start_loc)

% This function modifies matUINF to model a sinusoidal gust.

% 2017/06/04 - I-80W, Omaha, NE

% start_loc = repmat([-valGUSTSTART*valDELTIME*valUINF,0,0],size(fpg,1),1, size(fpg,3)); % Location (in meters) in global frame where gust starts

if flagGUSTMODE ~= 0
delx = start_loc - fpg; % Distance between DVE points and gust starting point

idx1 = delx(:,1) >= 0 & delx(:,1) <= valGUSTL;
idx2 = find(idx1 > 0);

idx3 = delx(:,1) >= 0 & delx(:,1) < 0.25*valGUSTL;
idx3_1 = find(idx3 > 0);
idx4 = delx(:,1) >= 0.25*valGUSTL & delx(:,1) < 0.75*valGUSTL;
idx4_1 = find(idx4 > 0);
idx5 = delx(:,1) >= 0.75*valGUSTL & delx(:,1) <= valGUSTL;
idx5_1 = find(idx5 > 0);

tau = delx(idx2,1)./valUINF;

gust_vel = zeros(size(fpg,1),1);

% Create gust velocity for sine wave gust
if flagGUSTMODE == 1
    
    if any(idx1) > 0
        tau = delx(idx2,1)./valUINF;
        gust_vel(idx2) = valGUSTAMP*(sin((2*pi*tau/(valGUSTL/valUINF))));
        matUINF(idx2,3) = matUINF(idx2,3) + (gust_vel(idx2));
    end

% Create gust velocity for 1-cosine gust
elseif flagGUSTMODE == 2
    
    if any(idx1) > 0
        gust_vel(idx2) = 0.5*valGUSTAMP*(1 - cos((2*pi*tau/(valGUSTL/valUINF))));
        matUINF(idx2,3) = matUINF(idx2,3) + (gust_vel(idx2));
%         gust_vel_old = gust_vel;
    end
    
% Create gust velocity for sharp edge gust
elseif flagGUSTMODE == 3
    
    if any(idx1) > 0
        gust_vel = valGUSTAMP;
        matUINF(idx2,3) = matUINF(idx2,3) + (gust_vel);
        gust_vel_old(idx2) = gust_vel;
    end
else
    
    disp('No gust mode exists for entered value')
    
end
    
else
    
    gust_vel = 0;
    
end


end
