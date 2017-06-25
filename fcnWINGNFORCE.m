function [valCL, valCLF, valCLI, valCY, valCYF, valCYI, valCDI, valE]= fcnWINGNFORCE(liftfree, liftind, sidefree, sideind, inddrag, matUINF, vecAREA, vecSPAN, vecSYM, valBETA, vecDVEVEHICLE, vecDVEWING)
%% Wing Normal Force
% this routine adds up the DVE's normal forces in order to compute the
% total wing normal forces/density and coefficients based on free stream

% INPUT:

% OUTPUT:
% valCL - Lift coefficient (freestream + induced)
% valCLF - lift coefficient due to freestream
% valCLI - Induced lift coefficient
% valCY - Side force coefficient (freestream + induced)
% valCYF - side force coefficient due to freestream
% valCYI - Induced side force coefficient

vehicles = unique(vecDVEVEHICLE(vecDVEWING > 0));
len = length(vehicles);

q = zeros(len,1);
valCL = zeros(1,len);
valCLF = zeros(1,len);
valCY = zeros(1,len);
valCLI = zeros(1,len);
valCDI = zeros(1,len);
AR = zeros(1,len);
valE = zeros(1,len);
valCYF = zeros(1,len);
valCYI = zeros(1,len);

for i = 1:len
    
    idxveh = vecDVEVEHICLE == i;
    %     q = 0.5.*sum(abs(vecUINF).^2,2).*valAREA;
    q(i) = 0.5.*mean((sqrt(sum(matUINF(idxveh,:).^2,2)).^2)).*vecAREA(i);
    
    ntfree = zeros(2,1);
    ntind = zeros(2,1);
    inddragsum = 0;
    
    %sum values from all elements
    ntfree(1) = sum(liftfree(idxveh));
    ntfree(2) = sum(sidefree(idxveh));
    
    ntind(1) = sum(liftind(idxveh));
    ntind(2) = sum(sideind(idxveh));
    
    inddragsum = sum(inddrag(idxveh));
    %double the force if we are using symmetry. This only works with sym for
    %the whole system
    if any(vecSYM) == 1 && valBETA ==0 %not sure why beta has to be zero
        disp('Symmetry is not currently working here! fcnWINGNFORCE');
        ntfree(1) = ntfree(1)*2; %why dont we double the side force?
        ntind(1) = ntind(1)*2;
        ntind(2) = ntind(2)*2;
        inddragsum = inddragsum*2;
    end
    
    %non-dimensionalize
    valCL(i) = (ntfree(1) + ntind(1))/q(i);
    valCLF(i) = ntfree(1)/q(i);
    
    valCY(i) = (ntfree(2) + ntind(2))/q(i);
    valCYF(i) = ntfree(2)/q(i);
    
    valCLI(i) = ntind(1)/q(i);
    valCYI(i) = ntind(2)/q(i);
    
    valCDI(i) = inddragsum/q(i);
    
    AR(i) = (vecSPAN(i)*vecSPAN(i))/vecAREA(i);
    valE(i) = (valCL(i)*valCL(i))/ (pi*AR(i)*valCDI(i));
end

end