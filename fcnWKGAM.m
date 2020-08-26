function vecWKGAM = fcnWKGAM(flagSTEADY, vecWKGAM, matCOEFF, vecDVETE, eta, valWNELE, valWSIZE)

if flagSTEADY == 1
    vecWKGAM = repmat([matCOEFF(vecDVETE>0,1) + ((eta.^2)./3).*matCOEFF(vecDVETE>0,3)], valWNELE/valWSIZE, 1);  
else
    vecWKGAM((end - valWSIZE + 1):end,1) = [matCOEFF(vecDVETE>0,1) + ((eta.^2)./3).*matCOEFF(vecDVETE>0,3)];
end

end