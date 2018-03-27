function [w_surf] = fcnSDVEVEL(fpg, SURF, INPU, FLAG)
%finds the velocity induced at a point due to all surface elements


% INPUT:
%  
% OUTPUT:
%   w_surf - [number points x 3] surface induced velocities

% List of sDVEs we are influencing from (each one for each of the fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:SURF.valNELE,len,1),[],1);

fpg = repmat(fpg,SURF.valNELE,1);

% DVE type 0 is a regular surface DVE
dvetype = zeros(length(dvenum),1);

[w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, SURF.matDVE, SURF.matVLST, SURF.matCOEFF, SURF.vecK, SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, SURF.vecDVELESWP, SURF.vecDVETESWP, INPU.vecSYM, FLAG.GPU);

w_surf = reshape(sum(reshape(w_ind', len*3, [])',1),3,[])';
end