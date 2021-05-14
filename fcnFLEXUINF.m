function [matUINF] = fcnFLEXUINF(matCENTER_t, matCENTER, valDELTIME, valTIMESTEP)

% lambda = 0.5;

matUINF = ((matCENTER_t(:,:,valTIMESTEP) - matCENTER)./(valDELTIME));
% matUINF = -((3*matCENTER - 4*matCENTER_t(:,:,valTIMESTEP) + matCENTER_t(:,:,valTIMESTEP-1))./(2*valDELTIME));

% yvel = (-matCENTER_t(:,2,valTIMESTEP-1) + 4*matCENTER_t(:,2,valTIMESTEP) - 3*matCENTER(:,2))./(2*valDELTIME);
% zvel = (-matCENTER_t(:,3,valTIMESTEP-1) + 4*matCENTER_t(:,3,valTIMESTEP) - 3*matCENTER(:,3))./(2*valDELTIME);

% matUINF(:,2) = yvel;
% matUINF(:,3) = zvel;
% matUINF = [(matCENTER_old(:,1) - matCENTER(:,1))./valDELTIME,(matCENTER(:,2) - matCENTER_old(:,2))./valDELTIME,(matCENTER(:,3) - matCENTER_old(:,3))./valDELTIME,];

% matUINF = (-matCENTER_old(:,:,valTIMESTEP-2) + 4*matCENTER_old(:,:,valTIMESTEP-1) - 3*matCENTER)./(2*valDELTIME);

end