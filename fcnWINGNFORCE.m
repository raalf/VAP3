function [valCL, valCLF, valCLI, valCY, valCYF, valCYI]= fcnWINGNFORCE(nfree,nind,liftfree,liftind,sidefree,sideind,vecUINF,valAREA,vecSYM,valBETA)
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


%  q=0.5*info.Uinf*info.Uinf*info.S
q = 0.5* sqrt(sum(abs(vecUINF).^2,2)) * valAREA;

ntfree(1) = 0;
ntfree(2) = 0;
ntind(1) = 0;
ntind(2) = 0;

%sum values from all elements
ntfree(1) = sum(liftfree);
ntfree(2) = sum(sidefree);

ntind(1) = sum(liftind);
ntind(2) = sum(sideind);

%double the force if we are using symmetry. This only works with sym for
%the whole system
if any(vecSYM) == 1 && valBETA ==0 %not sure why beta has to be zero 
    ntfree(1) = ntfree(1)*2; %why dont we double the ind lift?
    ntind(1) = ntfree(1)*2;
    ntind(2) = ntfree(2)*2;
end

%non-dimensionalize
valCL = (ntfree(1) + ntind(1))/q;
valCLF = ntfree(1)/q;

valCY = (ntfree(2) + ntind(2))/q;
valCYF = ntfree(2)/q;

valCLI = ntind(1)/q;
valCYI = ntind(2)/q;