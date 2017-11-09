function [a, b, c] = fcnDVEIND(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, flagGPU)

len = size(fpg,1); 

if flagGPU ~=1 || len < 6e+5
    [a, b, c] = fcnDVEIND_CHUNKS(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP);
else
    [a, b, c] = fcnDVEIND_GPU(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP);
end

end