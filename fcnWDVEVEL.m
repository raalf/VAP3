function [w_wake] = fcnWDVEVEL(fpg, valTIMESTEP, WAKE, SURF, FLAG)

% WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD, WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, ...
% WAKE.vecWDVEMCSWP, WAKE.vecWDVETESWP, WAKE.matWDVENORM, WAKE.matWVLST, WAKE.matWDVE, WAKE.valWNELE, ...
% WAKE.matWCENTER, WAKE.matWCOEFF, WAKE.vecWK

% List of wDVEs we are influencing from (each one for each of the fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:WAKE.valWNELE,len,1),[],1);

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

[w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, WAKE.matWDVE, WAKE.matWVLST, WAKE.matWCOEFF, WAKE.vecWK, WAKE.vecWDVEHVSPN, WAKE.vecWDVEHVCRD, WAKE.vecWDVEROLL, WAKE.vecWDVEPITCH, WAKE.vecWDVEYAW, WAKE.vecWDVELESWP, WAKE.vecWDVETESWP, SURF.vecDVESYM, FLAG.GPU);

w_wake = reshape(sum(reshape(w_ind', len*3, [])',1),3,[])';
end

