function [vecCLPDIST, vecCDPDIST, THRUSTDIST, TORQUEDIST] = fcnRVISCOUS3(flagVERBOSE, ...
    vecRPM,  valKINV, vecDVEHVCRD, vecN, vecM, ...
    vecDVELE, vecDVEPANEL, vecAIRFOIL, vecDISNORM, vecDVEAREA,...
    vecDVEROTOR, matUINF, matVLST, matDVE, matWUINF)
% This function applies a viscous correction using look up tables.
% OUTPUT
%   valCT - Viscous corrected thrust coeff
%   valCQ - Viscous corrected torque coeff
%   valCP - Viscous corrected power coeff

% Note: CN is normal to both the freestream and the spanwise direction
% Calculate velocity seen by section
% OLD velocity calculation
% tempXVEL = matUINF(:,1).*cos(vecTHETA)+matUINF(:,2).*sin(vecTHETA);
% vecV1 = sqrt(tempXVEL.^2 + matUINF(:,3).^2);

% Calculate chordline direction at midspan of each dve
avgle = (matVLST(matDVE(:,1),:)+matVLST(matDVE(:,2),:))./2;
avgte = (matVLST(matDVE(:,3),:)+matVLST(matDVE(:,4),:))./2;
tempDIF = avgte - avgle;
matCRDLINE = (tempDIF)./(repmat(sqrt(sum(tempDIF.^2,2)),[1,3]));
% Calculate velocity 
vecV = dot(matUINF + matWUINF, matCRDLINE,2);

% Calculate effective angle of attack
tempUINF = matUINF + matWUINF;
tempVELMAG = sqrt(tempUINF(:,1).^2 + tempUINF(:,2).^2 + tempUINF(:,3).^2);
tempCRDMAG = sqrt(matCRDLINE(:,1).^2 + matCRDLINE(:,2).^2 + matCRDLINE(:,3).^2);
%vecALPHAEFF = acos(dot(tempUINF, matCRDLINE,2)./(tempVELMAG.*tempCRDMAG));

[ledves, ~, ~] = find(vecDVELE > 0);
lepanels = vecDVEPANEL(ledves);

len = size(vecDISNORM,1);
tempDVEROTOR = ones(len,1);
idxdve = ledves(tempDVEROTOR(ledves) == 1);
idxpanel = lepanels(tempDVEROTOR(ledves) == 1);

m = vecM(idxpanel);

% Matrix of how much to add to an index to get the next chordwise element
tempm = repmat(vecN(idxpanel),1, m(1)).*repmat([0:m(1)-1],length(idxpanel),1);

% Which row index to add
rows = repmat(idxdve,1,m(1)) + tempm;

% Average velocities across chord
vecV = sum(vecV(rows),2)/(size(rows,2));

% CN = 2*(N/rho)/(V^2*S)
vecCNDIST = (sum(vecDISNORM(rows),2).*2)./(sum(vecDVEAREA(rows),2).*(vecV.^2));
vecALPHAEFF = vecCNDIST/(2*pi);
vecREDIST = vecV.*2.*sum(vecDVEHVCRD(rows),2)./valKINV;

% Different temp
vecCNDIST0 = vecCNDIST;
vecCDPDIST = zeros(size(rows,1),1);
vecCLPDIST = zeros(size(rows,1),1);
len = 0;
for j = 1:length(idxpanel)
    pan = idxpanel(j);
    airfoil = dlmread(strcat('airfoils/airfoil',num2str(vecAIRFOIL(pan)),'.dat'),'', 1, 0);

    HiRe = airfoil(end,4);
    LoRe = airfoil(1,4);

    alpha = vecALPHAEFF(len + j);

    if vecREDIST(len + j) > HiRe
        if flagVERBOSE == 1
            fprintf('\nRe higher than airfoil Re data')
        end
        Re2 = airfoil(end,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 1);
        cl_max_alpha = temp_var(end);
    elseif vecREDIST(len + j) < LoRe
        if flagVERBOSE == 1
            fprintf('\nRe lower than airfoil Re data');
        end
        Re2 = airfoil(1,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 1);
        cl_max_alpha = temp_var(end);
    else
        re1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 4);
        re1 = re1(end);
        cl_max_alpha1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 1);
        cl_max_alpha1 = cl_max_alpha1(end);

        temp_var = airfoil(airfoil(:,4) > vecREDIST(len + j),4);
        re2 = temp_var(1);
        temp_var = airfoil(airfoil(:,4) == (temp_var(1)), 1);
        cl_max_alpha2 = temp_var(end);

        cl_max_alpha = interp1([re1 re2],[cl_max_alpha1, cl_max_alpha2], vecREDIST(len + j));
    end
    % correcting the section cl if we are above cl_max
    if radtodeg(alpha) > cl_max_alpha
        if flagVERBOSE == 1
            fprintf('\nBlade Stall on Section %d, alpha = %f Re = %0.0f', j, radtodeg(alpha), vecREDIST(len + j))
        end
        %vecCNDIST0(len+j) = 0.825*cl_max; % setting the stalled 2d cl     
        
        % Make apply stall model using empirical equations
        % cn = cd,90*(sin(alpha_eff))/(0.56+0.44sin(alpha_eff))
        % ct = cd,0*cos(alpha_eff)/2
        % cd = cn*sin(alpha_eff)+ct*cos(alpha_eff)
        % Note: cd,90 = 2
        % Find cd_0
        temp = scatteredInterpolant(airfoil(:,4), airfoil(:,1), airfoil(:,3),'nearest');
        cd_0 = temp(vecREDIST(len+j),0);
        %alpha_eff = vecCNDIST/(2*pi);
        %alpha_eff = asin((vecCNDIST(len+j)/2)*0.56/(1-0.44*((vecCNDIST(len+j)/2))));
        
        cn = 2*sin(abs(vecALPHAEFF(len+j)))/(0.56+0.44*sin(abs(vecALPHAEFF(len+j))));
        ct = cd_0*cos(abs(vecALPHAEFF(len+j)))/2;
        vecCNDIST0(len+j)  = cn*cos(abs(vecALPHAEFF(len+j))) - ct*sin(abs(vecALPHAEFF(len+j)));
        
        vecCDPDIST(len + j) = cn*sin(abs(vecALPHAEFF(len+j))) + ct*cos(abs(vecALPHAEFF(len+j)));
        vecCLPDIST(len + j) = cn*cos(abs(vecALPHAEFF(len+j))) - ct*sin(abs(vecALPHAEFF(len+j)));
    else
        warning off
        F = scatteredInterpolant(airfoil(:,4), airfoil(:,1), airfoil(:,3),'nearest');
        vecCDPDIST(len + j, 1) = F(vecREDIST(len + j), radtodeg(alpha));
    end
end
% % Calculate viscous drag distribution
% vecDPDIST = 0.5*(vecCDPDIST.*((vecV).^2).*(sum(vecDVEAREA(rows),2)));
% 
% % Apply an appropriate directions
% tempUINFX = matUINF(:,1) + matWUINF(:,1);
% tempUINFY = matUINF(:,2) + matWUINF(:,2);
% tempUINFZ = matUINF(:,3) + matWUINF(:,3);
% matUINFAVG = [sum(tempUINFX(rows),2)/(size(rows,2)) sum(tempUINFY(rows),2)/(size(rows,2)) sum(tempUINFZ(rows),2)/(size(rows,2))];
% tempDIR = [matUINFAVG(:,1).*cos(vecTHETA(rows(:,1))) matUINFAVG(:,2).*sin(vecTHETA(rows(:,1))) matUINFAVG(:,3)];
% tempDIR = tempDIR./(sqrt(sum(tempDIR.^2,2)));
% matDPDIST = vecDPDIST.*tempDIR;

% Resolve to viscous thrust and torque
% THRUSTDIST = matDPDIST(:,3);
% TORQUEDIST = (dot(matDPDIST,[abs(cos(vecTHETA(rows(:,1)))) abs(sin(vecTHETA(rows(:,1)))) zeros(size(matDPDIST,1),1)],2)).*vecQARM(vecDVELE==1);
% POWERDIST = TORQUEDIST*2.*pi.*(vecRPM(vecDVEROTOR)./60);

THRUSTDIST = 1;
TORQUEDIST = 1;

% vecCNDISTDIF = vecCNDIST0 - vecCNDIST;
% vecDELNORMDISTP = 0.5*vecCNDISTDIF.*vecV.^2.*sum(vecDVEAREA(rows),2);
% difthrustP = vecDELNORMDISTP.*(en(:,3));
% diffsideP = vecDELNORMDISTP.*(dot(es,en,2));
% diffaxialP = vecDELNORMDISTP.*(dot(ea,en,2));
% Calculate coefficient
% valCTP = (sum(THRUSTDIST))/(((valRPM/60)^2)*((valDIA)^4));
% valCQP = (sum(TORQUEDIST))/(((valRPM/60)^2)*((valDIA)^5));
% valCPP = (sum(POWERDIST))/((valRPM/60)^3*(valDIA^5));


% Add viscous forces to thrust and power
% valCT = valCT + valCTP;
% valCQ = valCQP + valCQ;
% valCP = valCPP + valCP;
% if valTIMESTEP == 256
%     save('ViscForces225')
% end
% if valTIMESTEP == 228
%     save('ViscForces229')
% end
% if valTIMESTEP == 232
%     save('ViscForces233')
% end
% if valTIMESTEP == 236
%     save('ViscForces237')
% end
end