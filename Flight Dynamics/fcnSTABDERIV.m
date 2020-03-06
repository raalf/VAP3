function [TRIM, VEHI] = fcnSTABDERIV(COND, INPU, SURF, TRIM, VEHI, valTIMESTEP)

% Find dimensional and non-dimensional stability derivatives (convert all
% derivatives to /rad)

u0 = COND.vecVEHVINF;
theta0 = 0;
m = COND.vecVEHWEIGHT/9.81;
QS = 0.5*COND.valDENSITY*u0*u0*INPU.vecAREA;
Iy = 3.494;

% TRIM.valCD0 = TRIM.valCD;

% Aspect ratio of wing and tail
ARw = (INPU.vecSPAN*INPU.vecSPAN/INPU.vecAREA);
tempAR1 = SURF.matVLST(SURF.matDVE(SURF.vecWINGTYPE == 2,1),2);
tempAR2 = SURF.matVLST(SURF.matDVE(SURF.vecWINGTYPE == 2,2),2);
ARt = (2*(tempAR2(end)-tempAR1(1)))^2/(2*sum(SURF.vecDVEAREA(SURF.vecWINGTYPE == 2)));

CLalpha_t = 2*pi*ARt/(ARt + (2*(ARt+4)/(ARt+2))); % Lift curve slope for tail according to finite wing correction from McCormick (pg. 116)
CLalpha_w = 2*pi*ARw/(ARw + (2*(ARw+4)/(ARw+2))); % Lift curve slope for main wing according to finite wing correction from McCormick (pg. 116)


% Non-dimensional derivatives (Nelson Ch. 3), Dimensional derivatives
% (Etkin Ch. 4)
TRIM.Cw0 = m*9.81/QS;
% 

Vh = (SURF.valTBOOM*2*sum(SURF.vecDVEAREA(SURF.vecWINGTYPE == 2)))/(INPU.vecAREA*INPU.vecCMAC);

d_eps_dalpha = 2*CLalpha_w/(pi*ARw); % Rate of change of downwash with alpha (Nelson pg. 125)

TRIM.Cxu = -2*TRIM.valCD;
TRIM.Cxalpha = TRIM.valCL*(1 - 2*TRIM.CLalpha/(pi*TRIM.valE*ARw));
TRIM.Cxadot = 0;
TRIM.Cxq = 0;

TRIM.Czu = -2*TRIM.valCL*0;
TRIM.Czalpha = -(TRIM.CLalpha + TRIM.valCD);

TRIM.Czadot = -2*CLalpha_t*1*Vh*d_eps_dalpha;
TRIM.Czq = -2*CLalpha_t*1*Vh;
TRIM.Czde = -CLalpha_t*TRIM.tau*1*((2*sum(SURF.vecDVEAREA(SURF.vecWINGTYPE == 2)))/INPU.vecAREA);

% TRIM.Cmadot = TRIM.Czadot*(SURF.valTBOOM/INPU.vecCMAC);
% TRIM.Cmq = TRIM.Czq*(SURF.valTBOOM/INPU.vecCMAC);
TRIM.Cmu = 0;
TRIM.Cmq = -2*1*CLalpha_t*Vh*SURF.valTBOOM/INPU.vecCMAC;
TRIM.Cmadot = -2*1*CLalpha_t*Vh*SURF.valTBOOM*d_eps_dalpha/INPU.vecCMAC;

TRIM.CDalpha = 2*TRIM.valCL*TRIM.CLalpha/(pi*ARw*TRIM.valE);

% -------------------------- Nelson ---------------------------------------
% TRIM.Xu = -2*TRIM.valCD0*QS/(m*u0);
% TRIM.Xw = -(TRIM.CDalpha - TRIM.valCL)/(m*u0);
% 
% TRIM.Zu = -2*TRIM.valCL*QS/(m*u0);
% TRIM.Zw = -(TRIM.CLalpha + TRIM.valCD0)*QS/(m*u0);
% 
% TRIM.Mu = 0;
% TRIM.Mw = TRIM.Cmalpha*QS*INPU.vecCMAC/(u0*Iy);
% TRIM.Mq = TRIM.Cmq*INPU.vecCMAC/(2*u0);
% TRIM.Mwdot = TRIM.Cmadot*INPU.vecCMAC*QS*INPU.vecCMAC/(2*u0*u0*Iy);

% --------------------------- Etkin ---------------------------------------
TRIM.Czadot = -2*CLalpha_t*1*Vh*d_eps_dalpha;
TRIM.Czq = -2*CLalpha_t*1*Vh;
TRIM.Czde = -CLalpha_t*TRIM.tau*1*((2*sum(SURF.vecDVEAREA(SURF.vecWINGTYPE == 2)))/INPU.vecAREA);

TRIM.Cmadot = TRIM.Czadot*(SURF.valTBOOM/INPU.vecCMAC);
TRIM.Cmq = TRIM.Czq*(SURF.valTBOOM/INPU.vecCMAC);
TRIM.Cmu = 0;

TRIM.Xu = COND.valDENSITY*u0*INPU.vecAREA*sin(theta0) + 0.5*COND.valDENSITY*u0*INPU.vecAREA*TRIM.Cxu;
TRIM.Zu = -COND.valDENSITY*u0*INPU.vecAREA*TRIM.Cw0*cos(theta0) + 0.5*COND.valDENSITY*u0*INPU.vecAREA*TRIM.Czu;
TRIM.Mu = 0.5*COND.valDENSITY*u0*INPU.vecCMAC*INPU.vecAREA*TRIM.Cmu;

TRIM.Xw = 0.5*COND.valDENSITY*u0*INPU.vecAREA*TRIM.Cxalpha;
TRIM.Zw = 0.5*COND.valDENSITY*u0*INPU.vecAREA*TRIM.Czalpha;
TRIM.Mw = 0.5*COND.valDENSITY*u0*INPU.vecCMAC*INPU.vecAREA*TRIM.Cmalpha;

TRIM.Xwdot = 0;
TRIM.Zwdot = 0.25*COND.valDENSITY*INPU.vecCMAC*INPU.vecAREA*TRIM.Czadot; % Etkin
TRIM.Zwdot = -1*COND.valDENSITY*u0*CLalpha_t*(2*sum(SURF.vecDVEAREA(SURF.vecWINGTYPE == 2)))*SURF.valTBOOM*d_eps_dalpha/(2*u0);
TRIM.Mwdot = 0.25*COND.valDENSITY*INPU.vecCMAC*INPU.vecCMAC*INPU.vecAREA*TRIM.Cmadot*0;

TRIM.Xq = 0;
TRIM.Zq = 0.25*COND.valDENSITY*u0*INPU.vecCMAC*INPU.vecAREA*TRIM.Czq;
TRIM.Mq = 0.25*COND.valDENSITY*u0*INPU.vecCMAC*INPU.vecCMAC*INPU.vecAREA*TRIM.Cmq;

% Etkin state matrix
TRIM.matA = [TRIM.Xu/m, TRIM.Xw/m, 0, -9.81*cos(theta0);...
    TRIM.Zu/m, TRIM.Zw/m, (TRIM.Zq + m*u0)/m, (-m*9.81*sin(theta0))/m;...
    (1/Iy)*(TRIM.Mu + (TRIM.Mwdot*TRIM.Zu)/m), (1/Iy)*(TRIM.Mw + (TRIM.Mwdot*TRIM.Zw)/m),...
    (1/Iy)*(TRIM.Mq + (TRIM.Mwdot*(TRIM.Zq + m*u0))/m), (-TRIM.Mwdot*m*9.81*sin(theta0))/(Iy*m);...
    0, 0, 1, 0];

% Nelson state matrix
% A = [TRIM.Xu, TRIM.Xw, 0, -9.81; TRIM.Zu, TRIM.Zw, u0, 0; TRIM.Mu + TRIM.Mwdot*TRIM.Zu, TRIM.Mw + TRIM.Mwdot*TRIM.Zw, TRIM.Mq + TRIM.Mwdot*u0, 0; 0, 0, 1, 0];

[TRIM.vecEIGVEC, sys_eig] = eig(TRIM.matA); % Find eigenvectors and eigenvalues of the system
TRIM.vecEIGVAL = diag(sys_eig);

for i = 1:size(TRIM.matA,1)
    x_t(:,i) = TRIM.vecEIGVEC(:,i)*exp(TRIM.vecEIGVAL(i).*(valTIMESTEP*COND.valDELTIME));
end

TRIM.perturb(:,valTIMESTEP) = sum(x_t,2);
VEHI.matVEHUVW(1) = TRIM.matVEHUVW(1) - TRIM.perturb(1,valTIMESTEP);
VEHI.matVEHUVW(3) = TRIM.matVEHUVW(3) + TRIM.perturb(2,valTIMESTEP);
