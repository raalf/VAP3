function [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnELEVSWEEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC)

% Computing equivalent zero-lift line for a given elevator deflection,
% assuming an elevator chord of 25%. See McCormick. 
cf_c = 0.25; % Flap chord percentage (1 = All moving tail)
thetaf = acos(2*cf_c - 1);
TRIM.tau = 1 - (thetaf - sin(thetaf))/pi; % Flap effectiveness

de_seq = [(SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))/TRIM.tau; (SURF.vecDVEPITCH(SURF.idxTAIL(1)) - deg2rad(COND.vecVEHALPHA))/TRIM.tau + deg2rad(2)]; % Sequence is current elevator deflection and 2 degrees increase

for k = 1:length(de_seq)
    
    SURF.vecELEVANGLE = de_seq(k);
    
    if k > 1
        deltaEPS = TRIM.tau*(de_seq(k)-de_seq(k-1)); % Amount to rotate tail zero-lift line in rad
    else
        deltaEPS = TRIM.tau*(de_seq(k)-SURF.vecELEVANGLE); % Amount to rotate tail zero-lift line in rad
    end
    
    [SURF] = fcnTAILTRIM(SURF, FLAG, COND, deltaEPS, 1);

    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnVAP_TIMESTEP(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    [OUTP, COND, INPU, FLAG, MISC, SURF, TRIM, VEHI, VISC, WAKE] = fcnRESETVEHI(FLAG, COND, VISC, INPU, TRIM, VEHI, WAKE, SURF, OUTP, MISC);
    
    Cmde(k,1) = OUTP.vecVEHCM;
    CLde(k,1) = OUTP.vecCL(end);
    
end

% Line of best fit to get slope
Cm_de = polyfit(de_seq,Cmde,1);
TRIM.Cmde = Cm_de(1);

CL_de = polyfit(de_seq,CLde,1);
TRIM.CLde = CL_de(1);

end