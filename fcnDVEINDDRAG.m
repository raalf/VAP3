
function [inddrag]=fcnDVEINDDRAG(matCOEFF,matDVE,matVLST,vecUINF,vecDVEHVSPN,vecDVETE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, ...
    valWSIZE, valTIMESTEP,vecSYM,vecDVEWING )
% Induced dve drag. Function finds induced drag values on each te element. Outputs are not
% non-dimensionalized to q.



%% preliminary stuff
%te elements
idte = (vecDVETE == 3);

%number of te elements
numte = sum(idte);

%ABC for TE elements only

A = matCOEFF(idte,1);
B = matCOEFF(idte,2);
C = matCOEFF(idte,3);

%% FINDING INDUCED POINTS
% 80% of halfspan of TE elements only
eta8 = vecDVEHVSPN(idte).*0.8;

%TE vectors of TE elements only
s =( matVLST(matDVE(idte,3),:) -matVLST(matDVE(idte,4),:) )  ./ repmat((vecDVEHVSPN(idte).*2),1,3); %?? why is the S vector non-dim. to the span?

% element TE edge midpoint of TE elements only
xte = (matVLST(matDVE(idte,3),:) + matVLST(matDVE(idte,4),:))/2;

% we need to do a lot of work on the induced points
% New method of moving the TE points of the index (induced DVE) points.
% 17 Oct 2014. Bill B
% This method moves the points in the freestream direction into the plane
% passing through the TE of the inducer having freestream direction
% as the normal
tepoints = zeros(numte,3,3);

% find 3 points along TE of all TE elements
%first layer is left side, second layer is middle, third layer is right
%side
tepoints(:,:,1) = (xte + s.*repmat(-eta8,1,3)); %left side
tepoints(:,:,2) = xte ; %middle
tepoints(:,:,3) = (xte + s.*repmat(eta8,1,3)); %right ride

%% WORKING ON INDUCERS
% newest_row = [((valWNELE-valWSIZE)+1):1:valWNELE]';
% if on same wing!
% we won't use fcnWDVEVEL because we have to change the induced point for each
% column of wake elements, and change dve type for current timestep to 1. 
% So we will call fcnDVEVEL directly. So we have to set up all points/inducers

%need to repmat te points to get each tepoint numte times (induced)
%this will account for the induction of all the wake elements in the
%current timestep, on all the te points.
tepoints = repmat(tepoints,[numte,1,1]);

%dvenum is inducer
%need to keep the inducers index the same as the induced points
newest_row = [((valWNELE-valWSIZE)+1):1:valWNELE]';
dvenum = newest_row(repmat(1:valWSIZE,valWSIZE,1),:);

% now actually moving the point:
% vector from TE of each TE element to each point in tepoints
% order is as follows:
% influence of first dve on (1:numte), then influence of second dve
% on (1:numte), etc. 

delx  = tepoints-repmat(xte(repmat(1:numte,numte,1),:),[1 1 3]);

%project into freestream direction
temps = dot(delx,repmat(vecUINF,[size(delx,1) 1 3]),2);
tempb = repmat(temps,1,3,1).* repmat(vecUINF,[size(delx,1) 1 3]); %should this be normalized Uinf?

% original te point - tempb should be new te point
newtepoint = tepoints - tempb;

%end if on same wing, else newtepoint = oldtepoint

%we have now accounted for all the current timestep of wake elements, now repmat to
%account for remaining wake rows
%fpg is all points to go into DVEVEL
fpg = repmat(newtepoint,[valTIMESTEP,1,1]);
dvenum = repmat(dvenum,[valTIMESTEP,1,1]); %incorrect inducers index

mult = [1:valTIMESTEP]'; %need to renumber old timestep rows
multnew = repmat(mult,[valWSIZE*valWSIZE,1,1]);
multnew = sort(multnew);
dvenum = dvenum - repmat(valWSIZE,size(dvenum,1),1).*(multnew-1);
dvenum = repmat(dvenum,[1 1 3]);%correct inducers index

% take second dimension, move to bottom. then take third dimension and move
% to bottom
fpg = reshape(permute(fpg,[1 3 2]),[],3);
dvenum = reshape(permute(dvenum,[1 3 2]),[],1);

%dve type
dvetype = ones(length(dvenum),1);

dvetype(ismember(dvenum, newest_row)) = 1;%FW has this as type 1, but should be 2?

%setting singfct for post-trailing edge row to 0
tempwk = vecWK(dvenum);
tempwk(ismember(dvenum, newest_row)) = 0;

% Oldest row of wake DVEs are semi-infinite
oldest_row = [1:valWSIZE]';

if valTIMESTEP == 1
    dvetype(ismember(dvenum, oldest_row)) = 1;
else
    dvetype(ismember(dvenum, oldest_row)) = 3;
end

%get all velocities %need to set singfct = 0 for le row of elements!!!
[w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, matWDVE, matWVLST, matWCOEFF, tempwk, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, zeros(size(vecWDVELESWP)), zeros(size(vecWDVETESWP)), vecSYM);

% undo reshape and permute
w_total = permute(reshape(w_ind,[],3,3),[1 3 2]);

%add up influence of every dve on each point (every WSIZE value corresponds
%to the ind. vel on one point
w_wake(:,:,1) = reshape(sum(reshape(w_total(:,:,1)', valWSIZE*3, [])',1),3,[])';
w_wake(:,:,2) = reshape(sum(reshape(w_total(:,:,2)', valWSIZE*3, [])',1),3,[])';
w_wake(:,:,3) = reshape(sum(reshape(w_total(:,:,3)', valWSIZE*3, [])',1),3,[])';

%w_wake is [num tedves x 3 x k]

%% INTEGRATION
% //Kutta-Joukowski at 3 edges (left, center, right)
tempa = cross(w_wake,repmat(s,[1 1 3]),2);

gamma(:,1) = A - B.*eta8 + C.*eta8.*eta8;
gamma(:,2) = A;
gamma(:,3) =  A + B.*eta8 + C.*eta8.*eta8;

tempr = tempa.*repmat(permute(gamma,[1 3 2]),[1,3,1]);
% gamma1  = A - B.*eta8 + C.*eta8.*eta8;
% R1 = tempa(:,:,1).*repmat(gamma(:,1),1,3);

% gammao  = A;
% Ro = tempa(:,:,2).*repmat(gamma(:,2),1,3);

% gamma2  = A + B.*eta8 + C.*eta8.*eta8;
% R2 = tempa(:,:,3).*repmat(gamma(:,3),1,3);

% simpsons rule:
R(:,:)  = (tempr(:,:,1)+4.*tempr(:,:,2)+tempr(:,:,3)).*repmat(eta8,[1 3])./3;	
% R(:,:)  = (R1(:,:)+4*Ro(:,:)+R2(:,:)).*repmat(eta8,1,3)./3;	

% plus overhanging parts:
R(:,:) = R(:,:)+((7.*tempr(:,:,1)-8.*tempr(:,:,2)+7.*tempr(:,:,3)).*repmat(vecDVEHVSPN(idte)-eta8,[1 3])./3);
% R(:,:) = R(:,:)+((7.*R1(:,:)-8.*Ro(:,:)+7.*R2(:,:)).*repmat((vecDVEHVSPN(idte)-eta8),1,3)./3);

% R(:,1)  = (R1(:,1)+4*Ro(:,1)+R2(:,1)).*eta8./3;			%//Rx
% R(:,2)  = (R1(:,2)+4*Ro(:,2)+R2(:,2)).*eta8./3;			%//Ry
% R(:,3)  = (R1(:,3)+4*Ro(:,3)+R2(:,3)).*eta8./3;			%//Rz

% 		//plus overhanging parts
% R(:,1) = R(:,1)+((7.*R1(:,1)-8.*Ro(:,1)+7.*R2(:,1)).*(vecDVEHVSPN(idte)-eta8)./3); %//Rx
% R(:,2) = R(:,2)+((7.*R1(:,2)-8.*Ro(:,2)+7.*R2(:,2)).*(vecDVEHVSPN(idte)-eta8)./3); %//Ry
% R(:,3) = R(:,3)+((7.*R1(:,3)-8.*Ro(:,3)+7.*R2(:,3)).*(vecDVEHVSPN(idte)-eta8)./3); %//Rz
%% FORCES
inddrag(:,1) = dot(R,repmat(vecUINF,size(R,1),1),2);

end %end function