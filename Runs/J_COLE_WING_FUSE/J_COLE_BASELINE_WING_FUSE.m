clc
clear
warning off

seqALPHA = 1:1:12;

load('X57_FUSE.mat');
matFUSEORIG = [-3.7 0 -2.25];

% Running
clearvars -except seqALPHA vecCOLLECTIVE vecVEHVINF matFVLST matFDVE matFUSEORIG
filename = 'inputs/J_COLE_BASELINE_WING.vap';
for i = 1:length(seqALPHA)
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = seqALPHA(i);
    VAP_IN.valSTARTFORCES = 1;
    VAP_IN.valMAXTIME = 20;
    VAP_IN.RELAX = 0
    
    VAP_IN.matFDVE = matFDVE;
    VAP_IN.matFVLST = matFVLST;
    VAP_IN.matFUSEORIG = matFUSEORIG;
    VAP_IN.vecFDVEVEHICLE = ones(size(matFDVE,1),1);
    VAP_IN.vecFUSEVEHICLE = 1;
        
    OUTP(i) = fcnVAP_MAIN(filename, VAP_IN);
end

% save('VAP31_STEADY_VISCOUS_FIXED_SYM.mat')
