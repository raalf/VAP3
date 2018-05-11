function [out, ITER, ITEROUTP] = fcnBASELINE_OBJ()
% Define flight speed and conditions
vinf = 77.2;
rho = 1.007;
weightN = 13344.6648; % 3000 lb in N

N_prop = 1;

cd '../../'
seqALPHA = [2:10];
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valDELTIME = .25/vinf;
    VAP_IN.valSTARTFORCES = 30;
    VAP_IN.valMAXTIME = 30;
    VAP_IN.valSTARTFORCES = 3
    VAP_IN.valMAXTIME = 3
    WING_SWEEP(i) = fcnVAP_MAIN('X57_BASELINE_WING.vap', VAP_IN);
end
cd 'Runs/J_COLE_OPTIMIZATION/'

% Calculate CL at VINF and S&L flight
S    = WING_SWEEP(1).valAREA; % ref. platform area
CL   = weightN./(0.5*rho*vinf.^2*S);

% interpolate alpha to maintain steady level flight at VINF
% using wing only data
ALPHA = interp1([WING_SWEEP.vecCLv],[WING_SWEEP.vecVEHALPHA],CL);

% LD = interp1(borer(:,1),borer(:,2),vinf*1.94384); % get L/D from Borer Data
LD = 14;
CD = CL./(LD); % Calculate CD with Borer L/D Data
D  = 0.5*rho*vinf.^2.*CD*S; % Calulate drag force in Newton
thrust  = D/(2*N_prop); % Calculate Thrust force required from EACH PROP

%% Propeller Collective Sweep
cd '../../'
vecCOLLECTIVE = [-5:2:5];
for i = 1:length(vecCOLLECTIVE)
    VAP_IN = [];
    VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.valSTARTFORCES = 100;
    VAP_IN.valMAXTIME = 100;
    VAP_IN.valSTARTFORCES = 3
    VAP_IN.valMAXTIME = 3
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    PROP_SWEEP(i) = fcnVAP_MAIN('X57_BASELINE_PROP.vap', VAP_IN);
    %     view([90 90]);
end
cd 'Runs/J_COLE_OPTIMIZATION/'

%% Trim
ROTDIAM  = PROP_SWEEP(1).vecROTDIAM(1); % should be 1.524
ROTORRPM = PROP_SWEEP(1).vecROTORRPM(1); %should be 2250
rps      = ROTORRPM/60;

% Calculate CT required to maintain S&L flight
CT = thrust/(rho*rps^2*ROTDIAM^4);

% use scatteredInterpolant to avoid meshgrid
propCT   = [PROP_SWEEP.vecCT_AVG]';
propColl = [PROP_SWEEP.vecCOLLECTIVE]';
propVINF = repmat(vinf,size(propCT));

% meshgrids of 3d vinf,ct,collective data
propVINF = reshape(propVINF,[],length(unique(propVINF)));
propCT   = reshape(propCT,[],length(unique(propVINF)));
propColl = reshape(propColl,[],length(unique(propVINF)));

[~,maxCTidx] = max(propCT);
propCT(maxCTidx+1:end) = nan;

% idx = propColl<=0; % quick way to hack off the partially stalled propeller
% idx = ~isnan(propCT);
% F = scatteredInterpolant(propVINF(idx), propCT(idx), propColl(idx),'linear','nearest');

vecCOLLECTIVE = interp1(propCT, propColl, CT, 'linear', 'extrap');

vecCOLLECTIVE = vecCOLLECTIVE(~isnan(vecCOLLECTIVE));
CT = CT(~isnan(vecCOLLECTIVE));
vinf = vinf(~isnan(vecCOLLECTIVE));

ITER.maxIter = 5;
ITER.numCase = length(vinf);

ITER.Iteration = repmat((1:ITER.maxIter)',1,ITER.numCase);
ITER.CL        = nan(ITER.maxIter, ITER.numCase);
ITER.CT        = nan(ITER.maxIter, ITER.numCase);
ITER.CD        = nan(ITER.maxIter, ITER.numCase);
ITER.AOA       = nan(ITER.maxIter, ITER.numCase);
ITER.CLTV      = nan(ITER.maxIter, ITER.numCase);

seqALPHA = ALPHA;

for n = 1:ITER.maxIter
    
    if n == 2
        dCL1 = CL - ITER.CL(1,:);
        dCT1 = CT - ITER.CT(1,:);
        % New sets of AOA input for 2nd iteration in order to hit the targeted CL
        seqALPHA = interp1([WING_SWEEP.vecCLv],[WING_SWEEP.vecVEHALPHA],CL + dCL1, 'linear', 'extrap');
        % New sets of collective pitch input for 2nd iteration in order to hit the targeted CT
        vecCOLLECTIVE = interp1(propCT, propColl, CT + dCT1, 'linear', 'extrap');
    elseif n > 2
        % New sets of AOA input for the next iteration in order to hit the targeted CL
        seqALPHA = interp1(ITER.CL(1:n-1),ITER.AOA(1:n-1),CL,'linear','extrap');
        
        % New sets of collective pitch input input for the next iteration in order to hit the targeted CT
        vecCOLLECTIVE = interp1(ITER.CT(1:n-1),ITER.CLTV(1:n-1),CT,'linear','extrap');
        
    end
    
    ITER.AOA(n)  = seqALPHA;
    ITER.CLTV(n) = vecCOLLECTIVE;
    
    cd '../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA;
    VAP_IN.vecCOLLECTIVE = repmat(vecCOLLECTIVE, N_prop, 1);
    VAP_IN.vecVEHVINF = vinf;
    VAP_IN.valMAXTIME = 160;
    VAP_IN.valSTARTFORCES = VAP_IN.valMAXTIME-20;
    VAP_IN.valMAXTIME = 4
    VAP_IN.valSTARTFORCES = 2
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
    cd 'Runs/J_COLE_OPTIMIZATION/'
    
    % Write results
    ITER.CL(n,:) = OUTP.vecCL_AVG;
    ITER.CT(n,:) = nanmean(OUTP.vecCT_AVG);
    
    CDtemp = [OUTP.vecCD];
    CDtemp(isnan([OUTP.vecCLv])) = nan;
    ITER.CD(n,:) = nanmean(CDtemp,1);
    
    ITEROUTP(n).OUTP = OUTP;
    
end

TRIMMED = false;

dCT = abs((ITER.CT(end) - CT)./CT);
dCL = abs((ITER.CL(end) - CL)./CL);

if dCT <= 0.01 && dCL <= 0.01
    TRIMMED = true;
end


%% ANALYZE RESULTS
if TRIMMED == true
    out = sum(2.*OUTP.vecCP_AVG).*((2250/60).^3).*(1.524.^5).*OUTP.valDENSITY;
    
else
    out = 1e+20;
end

end

