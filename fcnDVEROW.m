function [matROWS] = fcnDVEROW(ledves, SURF, INPU)

% This functions determines what DVEs are in which row or column across the
% span of the wing. matROWS will be a matrix that is NELE x INPU.vecM.

lepanels = SURF.vecDVEPANEL(ledves);

% Determine DVEs in each spanwise station
for i = 1:max(SURF.vecDVEWING)

	idxdve = ledves(SURF.vecDVEWING(ledves) == i);
	idxpanel = lepanels(SURF.vecDVEWING(ledves) == i);

    m = INPU.vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
    
    matROWS{i} = repmat(idxdve,1,m) + double(tempm);

end

end