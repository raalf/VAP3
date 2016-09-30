function [w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM)
% This function takes in DVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 

% T.D.K 2016-09-11 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

len = length(dvenum);

[a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM);

D = [a b c];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

w_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,3,[]),3,1,1),2),[2 1 3]);
w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);

end