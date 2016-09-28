function [vecUINF] = fcnUINFWING(valALPHA, valBETA)

% INPUT:
%   valALPHA - Angle of attack for this move (radians)
%   valBETA - Sideslip angle for this move (radians)
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
% OUTPUT:
%   vecUINF - vector of freestream velocity

uinf = 1;

vecUINF = [uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)];

end

