function [w_surf] = fcnSDVEVEL(fpg, valNELE, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM)
%finds the velocity induced at a point due to all surface elements


% INPUT:
%  
% OUTPUT:
%   w_wsurf - length(fpg) x 3  induced velocities

% List of sDVEs we are influencing from (each one for each of the fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);

fpg = repmat(fpg,valNELE,1);

% DVE type 0 is a regular surface DVE
dvetype = zeros(length(dvenum),1);

[w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM);

w_surf = reshape(sum(reshape(w_ind', len*3, [])',1),3,[])';
end