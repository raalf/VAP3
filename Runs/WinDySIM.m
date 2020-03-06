clc
clear
warning off

filename = 'inputs/WinDySIM.vap';

VAP_IN = [];

[FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n] = fcnVAPINIT(filename,VAP_IN);

COND.valTRIMSTEP = 1;

COND.valFULLTRIMSTEP = 0;

FLAG.FULLTRIM = 0; % Flag indicating whether the vehicle is fully trimmed with tail and aeroelastic deflections

% Initialize trim tolerances
trim_tol_def = 100;
trim_tol_twist = 100;

tic;

while FLAG.FULLTRIM == 0
    
    % Compute static aeroelastic response of wing
    while FLAG.STATICAERO == 0
        
        COND.valFULLTRIMSTEP = COND.valFULLTRIMSTEP + 1;

    	% Set flag to not trim tail
    	FLAG.TRIM = 0;
        
        [SURF, INPU, VEHI, OUTP, COND, FLAG] = fcnELASTICVEHI(SURF, INPU, VEHI, OUTP, COND, FLAG);

    	[OUTP, SURF, COND, INPU, MISC, VEHI, VISC, FLAG, matD, vecR] = fcnVAP_MAIN(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n);
        
        OUTP.matDEFGLOBTRIM(COND.valFULLTRIMSTEP,:) = OUTP.matDEFGLOB(end,:);
        OUTP.matTWISTGLOBTRIM(COND.valFULLTRIMSTEP,:) = OUTP.matTWISTGLOB(end,:);
        
    end

    % Set flag to trim tail once static aeroelastic convergence has been reached
    FLAG.TRIM = 1;
    FLAG.STIFFWING = 1;

	[OUTP, SURF, COND, INPU, MISC, VEHI, VISC, FLAG, matD, vecR] = fcnVAP_MAIN(FLAG, COND, VISC, INPU, VEHI, WAKE, SURF, OUTP, MISC, matD, vecR, n);
    
    % Compute change in aeroelastic response to determine complete
    % converged/trimmed solution
    if COND.valFULLTRIMSTEP > 1 && FLAG.TRIMMED == 1
        
        trim_tol_def = 100*(OUTP.matDEFGLOBTRIM(COND.valFULLTRIMSTEP,end) - OUTP.matDEFGLOBTRIM(COND.valFULLTRIMSTEP-1,end))/OUTP.matDEFGLOBTRIM(COND.valFULLTRIMSTEP-1,end);
        trim_tol_twist = 100*(OUTP.matTWISTGLOBTRIM(COND.valFULLTRIMSTEP,end) - OUTP.matTWISTGLOBTRIM(COND.valFULLTRIMSTEP-1,end))/OUTP.matTWISTGLOBTRIM(COND.valFULLTRIMSTEP-1,end);
        
        if trim_tol_def < 2 && trim_tol_twist < 2
            FLAG.FULLTRIM = 1;
        end
        
    end
    
    % Reset static aeroelastic convergence to compute new configuration
    % with new tail angle
    if FLAG.TRIMMED == 1 && FLAG.FULLTRIM == 0
        FLAG.STATICAERO = 0;
        FLAG.TRIMMED = 0;
        FLAG.STIFFWING = 0;
    end

    COND.valTRIMSTEP = COND.valTRIMSTEP + 1;
        
    if FLAG.FULLTRIM == 1
        tail_incidence = (SURF.vecDVEPITCH(SURF.vecDVEWING == 2)*180/pi-COND.vecVEHALPHA);
        fprintf('\nSuccess! Vehicle trimmed with a tail incidence of %.2f deg.\n', tail_incidence(1))
    end
    
end

toc;