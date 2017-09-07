function [sign] = etaSign(EdgeNumber)
% etaSign Gives a sign based on the Edge number, if the edge number
% is one, then it returns -1. If the edge number is 2, it returns 1. 

if EdgeNumber == 1
    sign = -1;
elseif EdgeNumber == 2
    sign = 1;
end

end

