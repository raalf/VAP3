function [w_wake, WAKE] = fcnWDVEVEL_SHORTCUT(influencee, valTIMESTEP, WAKE, SURF, FLAG)
% INPUT:
%   influencee - DVE number of the influencee (influence is on this DVE's control point
if valTIMESTEP == 1
    WAKE.matWINDUC_REF = nan(WAKE.valWSIZE, SURF.valNELE, 1);
elseif valTIMESTEP == 2
    WAKE.matWINDUC_REF = nan(WAKE.valWSIZE.*2, SURF.valNELE, 1);
else
    % Upper 2 rows (oldest wake elements) are nan, because these are semi-infinite sheets and have to be recalculated
    WAKE.matWINDUC_REF = [nan(WAKE.valWSIZE.*2, SURF.valNELE, 1); WAKE.matWINDUC_REF(WAKE.valWSIZE+1:end, :)];
end

% List of wDVEs we are influencing from (each one for each of the fieldpoints)
len_original = length(influencee);
len = length(influencee);
dvenum = reshape(repmat(1:WAKE.valWNELE,len,1),[],1);
influencee = repmat(influencee, WAKE.valWNELE, 1);

fpg = SURF.matCENTER;
fpg = repmat(fpg,WAKE.valWNELE,1);

if FLAG.STEADY == 1
    % DVE type 1 is a regular wake DVE
    dvetype = ones(length(dvenum),1);
    
    % Newest row of wake DVEs have a filament at the leading edge
    if FLAG.TRI == 1
        newest_row = sort([WAKE.valWNELE:-2:WAKE.valWNELE-WAKE.valWSIZE*2+1]-1)';
    else
        newest_row = [((WAKE.valWNELE-WAKE.valWSIZE)+1):1:WAKE.valWNELE]';
    end
    
    dvetype(ismember(dvenum, newest_row)) = 2;
    
    % Oldest row of wake DVEs are semi-infinite
    if FLAG.TRI == 1
        oldest_row = [2:2:WAKE.valWSIZE*2];
    else
        oldest_row = [1:WAKE.valWSIZE]';
    end
    
    if valTIMESTEP == 1 && FLAG.TRI == 0
        dvetype(ismember(dvenum, oldest_row)) = -3;
    elseif valTIMESTEP > 0 && FLAG.TRI == 0
        dvetype(ismember(dvenum, oldest_row)) = 3;
    elseif FLAG.TRI == 1 % last row of wake elements in a tri-wake never has a bound vortex
        dvetype(ismember(dvenum, oldest_row)) = 3;
    end
    
elseif FLAG.STEADY == 0 || 2
    % DVE type 1 is a regular wake DVE
    dvetype = zeros(length(dvenum),1);
    
    % Oldest row of wake DVEs are semi-infinite
    if FLAG.TRI == 1
        oldest_row = [2:2:WAKE.valWSIZE*2];
    else
        oldest_row = [1:WAKE.valWSIZE]';
    end
    
    dvetype(ismember(dvenum, oldest_row)) = -3;
    
end

%% Finding stored calculations
len = length(dvenum);
linearInd = sub2ind(size(WAKE.matWINDUC_REF), dvenum, influencee);
stored = ~isnan(WAKE.matWINDUC_REF(linearInd));
a = nan(len,3);
b = nan(len,3);
c = nan(len,3);

if valTIMESTEP > 1 && ~isempty(WAKE.matWINDUC)
    a(stored,:) =  WAKE.matWINDUC(WAKE.matWINDUC_REF(linearInd(stored)),:,1);
    b(stored,:) =  WAKE.matWINDUC(WAKE.matWINDUC_REF(linearInd(stored)),:,2);
    c(stored,:) =  WAKE.matWINDUC(WAKE.matWINDUC_REF(linearInd(stored)),:,3);
end

%% Finding unstored computations
[a(~stored,:), b(~stored,:), c(~stored,:)] = fcnDVEINF(dvenum(~stored), dvetype(~stored), fpg(~stored,:), WAKE.vecWK, ...
    WAKE.matWDVE, WAKE.matWVLST, WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD, WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, ...
    WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, WAKE.vecWDVETESWP, SURF.vecDVESYM, FLAG.GPU);

% Saving the influence coefficients
if valTIMESTEP > 1
    if isempty(WAKE.matWINDUC)
        idx = 1:len;
    else
        idx = size(WAKE.matWINDUC, 1) + 1:size(WAKE.matWINDUC,1) + sum(~stored);
    end
    WAKE.matWINDUC_REF(sub2ind(size(WAKE.matWINDUC_REF), dvenum(~stored), influencee(~stored))) = idx;
    WAKE.matWINDUC(idx,:,1) = a(~stored,:);
    WAKE.matWINDUC(idx,:,2) = b(~stored,:);
    WAKE.matWINDUC(idx,:,3) = c(~stored,:);
end

%% Finding the induced velocities
D = [a b c];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

w_ind = permute(sum(D.*repmat(reshape(WAKE.matWCOEFF(dvenum,:)',1,3,[]),3,1,1),2),[2 1 3]);
w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);
w_wake = reshape(sum(reshape(w_ind', len_original*3, [])',1),3,[])';

end

