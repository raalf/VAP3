function [vecCL] = fcnUNSTEADYWRAPPER(vecDVENFREE,vecUINF,vecDVEHVCRD,valAREA,gamma_old,...
    valMAXTIME,valGUSTSTART,valDELTIME,vecSYM,valBETA,en_t)
% The unsteady wrapper, first derived on Embraer E175, AC7665, 33000 ft.
% Implemented: Fishy AirBnB, Dallas, TX 75219

q = 0.5* sqrt(sum(abs(vecUINF).^2,2)) * sqrt(sum(abs(vecUINF).^2,2)) * valAREA;

gamma_old(:,valMAXTIME+1) = gamma_old(:,valMAXTIME); % Sketchy....standy by

% idx1 = vecDVELE == 1; %index of LE vectors (will be the same)
% 
% % find vector across element (should already have this somewhere...)
% % for first spanwise row, vector is LE vect, for all other spanwise rows,
% % vector is halfchord vect.
% 
% % vector of bound vortex along leading edge for LE row
% % tempa(idx1,1) = tan(vecDVELESWP(idx1));
% s = zeros(valNELE,3);
% s(idx1,:) =( matVLST(matDVE(idx1,2),:) -matVLST(matDVE(idx1,1),:) )  ./ repmat((vecDVEHVSPN(idx1).*2),1,3); %?? why is the S vector non-dim. to the span?
% 
% % UxS
% tempb = cross(matUINF,s);
% 
% % norm(UxS)
% uxs = sqrt(sum(abs(tempb).^2,2));
% 
% % eN = tempa.*(1/UxS);
% en = tempb.*repmat((1./uxs),1,3);

dGammadt(:,valGUSTSTART:valMAXTIME) = (gamma_old(:,valGUSTSTART+1:valMAXTIME+1) - gamma_old(:,valGUSTSTART:valMAXTIME))./(valDELTIME); % Forward difference up until last time step
% dGammadt(:,valGUSTSTART:valMAXTIME) = (gamma_old(:,valGUSTSTART+1:valMAXTIME+1) - gamma_old(:,valGUSTSTART-1:valMAXTIME-1))./(2*valDELTIME); % Central difference up until last time step
% dGammadt(:,valGUSTSTART:valMAXTIME) = (3*gamma_old(:,valGUSTSTART:valMAXTIME) - 4*gamma_old(:,valGUSTSTART-1:valMAXTIME-1) + gamma_old(:,valGUSTSTART-2:valMAXTIME-2))./(2*valDELTIME);
dGammadt(:,valMAXTIME) = (3*gamma_old(:,valMAXTIME) - 4*gamma_old(:,valMAXTIME-1) + gamma_old(:,valMAXTIME-2))./(2*valDELTIME);
dGammadt(:,valGUSTSTART) = (3*gamma_old(:,valGUSTSTART) - 4*gamma_old(:,valGUSTSTART-1) + gamma_old(:,valGUSTSTART-2))./(2*valDELTIME);

vecDVENFREE = vecDVENFREE + dGammadt;
vecDVENFREE = reshape(vecDVENFREE,[size(vecDVENFREE,1) 1 valMAXTIME]);
liftfree = vecDVENFREE.*en_t(:,3,:);
ntfree = sum(liftfree,1);

if any(vecSYM) == 1 && valBETA ==0 %not sure why beta has to be zero 
    ntfree = ntfree*2; %why dont we double the side force?
end

vecCL = (ntfree)./q;

vecCL = reshape(vecCL,[valMAXTIME,1,1]);

end