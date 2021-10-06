function [w_total] = fcnINDVEL(fpg, valTIMESTEP, SURF, WAKE, INPU, FLAG)

%returns the induced velocity at a point due to all surface and wake
%elements


% INPUT:
%   pass in a vector of points (num points x 3)
% OUTPUT:
%   w_total - numpoints x 3  induced velocities

% velocities from surface elements
[w_surf] = fcnSDVEVEL(fpg, SURF, INPU, FLAG);

% velocities from wake elements
w_wake = fcnWDVEVEL(fpg, valTIMESTEP, WAKE, SURF, FLAG);

% add
w_total = w_surf + w_wake;
% w_total = w_surf;

