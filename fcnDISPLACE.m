function [x] = fcnDISPLACE(vup, vnow, vdown, x, valDELTIME)

% This function takes a point and moves it based on
% the velocity at the locations around it

% It will NOT work for rotors now, for that, we will
% need to consider the effects of UINF 

u = 0.25.*(vdown + 2.*vnow + vup);

x = x + (u.*valDELTIME);

end

