function [OUTP] = fcnVEHENERGY(INPU, COND, SURF, OUTP, VEHI, FLAG, valTIMESTEP)
tot_mass = (sum(2*SURF.vecVEHMASS,1) + sum(VEHI.vecFUSEMASS,1) + VEHI.vecWINGMASS(2) + VEHI.vecPAYLMASS);
% Index of leading edge DVEs on beam
idx = SURF.vecDVELE(SURF.idxFLEX)==1;
idx1 = find(idx == 1);

% ele_vel = vecnorm((SURF.matCENTER(idx1,:) - SURF.matCENTER_t(idx1,:,valTIMESTEP))./COND.valDELTIME,2,2); % Calculate velocity of each DVE element
% ele_vel = ((OUTP.matDEFGLOB(valTIMESTEP,:) - OUTP.matDEFGLOB(valTIMESTEP-1,:))./COND.valDELTIME)';
% ele_vel = interp1(SURF.matBEAMLOC(:,2,valTIMESTEP),ele_vel,SURF.matCENTER(idx1,2));
fuse_vel = COND.vecVEHVINF; % Velocity of fuselage/tail system
% OUTP.ele_vel(:,valTIMESTEP) = ele_vel;

% OUTP.work(valTIMESTEP,1) = OUTP.GlobForce(valTIMESTEP,3)*(OUTP.vecCGLOC(valTIMESTEP,3)-OUTP.vecCGLOC(valTIMESTEP-1,3));

if FLAG.STIFFWING == 1
    ele_vel = vecnorm(repmat(OUTP.matGLOBUVW(valTIMESTEP,:),INPU.valNSELE,1),2,2);
else
    ele_vel = vecnorm([repmat(OUTP.matGLOBUVW(valTIMESTEP,1),length(OUTP.valUDOT(:,valTIMESTEP)),1), repmat(OUTP.matGLOBUVW(valTIMESTEP,2),length(OUTP.valUDOT(:,valTIMESTEP)),1),...
        OUTP.matGLOBUVW(valTIMESTEP,3) + OUTP.valUDOT(:,valTIMESTEP)],2,2);
end

% Kinetic energy of wing elements and fuselage/tail system
OUTP.ele_kinen_temp(:,valTIMESTEP) = 0.5.*[0; SURF.vecVEHMASS].*ele_vel.*ele_vel;
OUTP.ele_kinen(valTIMESTEP,1) = sum(OUTP.ele_kinen_temp(:,valTIMESTEP),1);
% OUTP.ele_kinen(valTIMESTEP,1) = trapz(SURF.matCENTER(idx1,2),OUTP.ele_kinen_temp(:,valTIMESTEP));
% OUTP.fuse_kinen(valTIMESTEP,1) = 0.5.*VEHI.vecFUSEMASS.*fuse_vel.*fuse_vel;

% Gravitational potential energy of wing elements and fuselage/tail system
OUTP.ele_poten_temp(:,valTIMESTEP) = SURF.vecVEHMASS.*9.81.*SURF.vecWINGCG(:,3);
OUTP.ele_poten(valTIMESTEP,1) = sum(OUTP.ele_poten_temp(:,valTIMESTEP),1);
% OUTP.fuse_poten(valTIMESTEP,1) = VEHI.vecFUSEMASS.*9.81.*VEHI.vecFUSECG(:,3);

%% Calculate curvature of deflected wing
Ixx = SURF.vecVEHMASS.*SURF.matCENTER(idx1,2).*SURF.matCENTER(idx1,2); % Mass moment of inertia about fuselage axis

% OUTP.ele_roten_temp(:,valTIMESTEP) = Ixx.*SURF.matBEAMOMEGA(:,1,valTIMESTEP).*SURF.matBEAMOMEGA(:,1,valTIMESTEP);
% OUTP.ele_roten(valTIMESTEP,1) = sum(OUTP.ele_roten_temp(:,valTIMESTEP),1);

% Perform second derivative of beam displacement using central finite
% differencing
curv = (OUTP.matDEF(end,4:end-1) - 2*OUTP.matDEF(end,3:end-2) + OUTP.matDEF(end,2:end-3))'./([INPU.valDY(1); INPU.valDY].*[INPU.valDY(1); INPU.valDY]);

M = -INPU.matEIx(:,1).*curv;

% V = (M(3:end) - M(1:end-2))./(2.*INPU.valDY(2:end));
% V = [V; (3*M(end) - 4*M(end-1) + M(end-2))./(2.*INPU.valDY(end))];
% V = [(-3*M(1) + 4*M(2) - M(3))./(2*INPU.valDY(1)); V];
OUTP.vecFUSEV(valTIMESTEP,1) = trapz(SURF.vecSTRUCTSPNDIST,OUTP.vecBEAMFORCE);
OUTP.vecFUSEM(valTIMESTEP,1) = trapz(SURF.vecSTRUCTSPNDIST,OUTP.vecBEAMMOM);
U = trapz(SURF.vecSTRUCTSPNDIST,(M.*M)./(2.*INPU.matEIx(:,1)));

U = trapz(SURF.vecSTRUCTSPNDIST,((-INPU.matEIx(:,1).*curv).^2)./(2.*INPU.matEIx(:,1)));

%% Wing element elastic energy
OUTP.temp_ele_U(:,valTIMESTEP) = 0.5.*INPU.matEIx(:,1).*curv.^2; % Determine integrand of elastic energy equation
% OUTP.temp_ele_U(:,valTIMESTEP) = (M.*M)./(2.*EI); % Determine integrand of elastic energy equation

% Integrate between adjacent elements to determine elastic energy of
% section
for i = 2:length(curv)
%     OUTP.ele_U_temp(i,valTIMESTEP) = trapz(SURF.matCENTER(idx1(i-1:i),2),OUTP.temp_ele_U(i-1:i,valTIMESTEP));
    OUTP.ele_U_temp(i,valTIMESTEP) = trapz(SURF.vecSTRUCTSPNDIST(i-1:i,1),OUTP.temp_ele_U(i-1:i,valTIMESTEP));
end

OUTP.ele_U(valTIMESTEP,1) = sum(OUTP.ele_U_temp(:,valTIMESTEP),1);

% OUTP.work(valTIMESTEP,1) = trapz(SURF.vecSTRUCTSPNDIST,OUTP.vecBEAMFORCE.*(OUTP.matDEFGLOB(valTIMESTEP,:)-OUTP.matDEFGLOB(valTIMESTEP-1,:))');

%% Compute total energy (element approach)
OUTP.tot_en(valTIMESTEP,1) = 2*OUTP.ele_U(valTIMESTEP,1) + 2*OUTP.ele_poten(valTIMESTEP,1) + 2*OUTP.ele_kinen(valTIMESTEP,1); % Total energy

OUTP.vecZE(valTIMESTEP,1) = OUTP.tot_en(valTIMESTEP,1)./(tot_mass*9.81); % Energy altitude

%% Compute total energy (system approach)
% OUTP.U(valTIMESTEP,1) = trapz(SURF.matCENTER(idx1,2),OUTP.temp_ele_U(:,valTIMESTEP)); % Elastic energy of wing
OUTP.U(valTIMESTEP,1) = trapz(SURF.vecSTRUCTSPNDIST,OUTP.temp_ele_U(:,valTIMESTEP)); % Elastic energy of wing
OUTP.poten(valTIMESTEP,1) = tot_mass*9.81*INPU.vecVEHCG(:,3);
OUTP.kinen(valTIMESTEP,1) = 0.5.*tot_mass.*COND.vecVEHVINF.*COND.vecVEHVINF;
OUTP.tot_en_1(valTIMESTEP,1) = OUTP.kinen(valTIMESTEP,1) + OUTP.poten(valTIMESTEP,1) + 2*OUTP.U(valTIMESTEP,1);
OUTP.vecZE(valTIMESTEP,2) = OUTP.tot_en_1(valTIMESTEP,1)./(tot_mass*9.81);