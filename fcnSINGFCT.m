function [vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN)
% Each DVE is assigned a singularity factor, which is 1% of the half-span
% of the wingtip DVE. That means DVEs on different wings have different
% singularity factors. This factor is used to help with log singularities in
% vortex sheet induction.

% INPUT:
%   valNELE - total number of DVEs
%   vecDVEWING - valNELE x 1 vector of wing number for each DVE
%   vecDVETIP - valNELE x 1 vector with local edges where a DVE is at the wingtip
%   vecDVEHVPSN - valNELE x 1 vector of DVE half spans
% OUTPUT:
%   vecK - valNELE x 1 vector of singularity factors

[~, ia, ~] = unique([vecDVEWING vecDVETIP],'rows');

idx = vecDVETIP(ia) > 0;

k = vecDVEHVSPN(ia(idx));
wings = vecDVEWING(ia(idx));

% To apply singularity factor, I create a matrix of a bunch of 0's and 1's, where
% each DVE has a 1 in a column. The column that 1 is in corresponds with the wing 
% number, so by repmatting the k for each wing, and multiplying them together,
% I get the proper K value for each DVE based on what wing it belongs to.
% Then it is sorted in each row to get rid of the zeros
vecK = sort([repmat(vecDVEWING, 1, max(wings)) == repmat(wings', valNELE,1)].*repmat(k', valNELE,1),2,'descend');

vecK = 0.01.*vecK(:,1);

end

