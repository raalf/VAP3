function [w_ind] = fcnDVEVEL_OL(dvenum, fpg, dvetype, matDVE, matVLST, matCOEFF, vecK, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM, matESHEETS)
% This function takes in DVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 

% T.D.K 2016-09-11 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

len = length(dvenum);

[a, b, c, d, e, f] = fcnDVEINF_OL(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM, matESHEETS);

D = [a b c d e f];
% D = [d e f a b c];
D = reshape(reshape(D', 1, 18, []), 3, 6, len);

w_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);
w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);

end