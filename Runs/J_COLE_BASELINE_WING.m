clc
clear
warning off

seqALPHA = 10;
vecCOLLECTIVE = seqALPHA.*0;

filename = 'inputs/J_COLE_BASELINE_WING.vap';
for i = 1:length(vecCOLLECTIVE)
    OUTP(i) = fcnVAP_MAIN(filename, seqALPHA(i), vecCOLLECTIVE(i));
end
