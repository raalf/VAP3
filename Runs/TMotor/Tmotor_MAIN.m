clc
clear
warning off

filename = 'inputs/TMotor.vap';
vecJ_1 = 0.1:0.1:1;
vecALPHA = -1*[-90 -30 -15 -5 0 5 15 30 60 90];
vecVEHFPA = -1*vecALPHA;
vecJ_2 = 0.1:0.05:0.5;
% VAP_IN = [];
% OUTP = fcnVAP_MAIN(filename, VAP_IN);

for i = 1:length(vecALPHA)
    if vecALPHA(i) <= -60
        seqJ = vecJ_2;
    else 
        seqJ = vecJ_1;
    end
    fpa = vecVEHFPA(i); 
    alpha = vecALPHA(i);
    vecALPHA(i)
    parfor j = 1:length(seqJ)
        t = datestr(now)
        fprintf('\n[%s]\t\t Running: Alpha = %0.1f, J = %0.3f\n',t, fpa, seqJ(j))
        VAP_IN = [];
%         VAP_IN.valSTARTFORCES = 139;
        VAP_IN.vecCOLLECTIVE = 2.6;
        
        VAP_IN.vecROTORRPM = 3000;
        VAP_IN.RELAX = 1;
        
        VAP_IN.vecVEHFPA = fpa;
        VAP_IN.vecVEHALPHA = alpha;
        VAP_IN.vecVEHVINF = seqJ(j)*0.4572*3000/60;
        
        OUTP(i,j) = fcnVAP_MAIN(filename, VAP_IN);
        t = datestr(now)
        fprintf('\n[%s]\t\t Complete: Alpha = %0.1f, J = %0.3f\n',t, fpa, seqJ(j))
    end
end

save('Runs\TMotor\VAP3_3000RPM_Relaxed')