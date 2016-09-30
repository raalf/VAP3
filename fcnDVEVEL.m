function [w_ind] = fcnDVEVEL(dvenum, fpg, matDVE, matDVENORM, matVLST)
% This function takes in DVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 

% T.D.K 2016-09-11 6075 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

len = length(dvenum);

[a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM)

D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);


q_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);

q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

end