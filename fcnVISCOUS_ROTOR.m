function [matCDPDIST] = fcnVISCOUS_ROTOR(flagVERBOSE, valKINV, ...
    vecDVEHVCRD, vecN, vecM, vecDVELE, vecDVEPANEL, cellAIRFOIL, vecDISNORM, vecDVEAREA, matUINF, matVLST, matDVE, matWUINF)
% This function applies a viscous correction to rotors using look up tables
% and applies a high angle stall model.
% 
% OUTPUT
%   matCDPDIST - CDP with direction accounted for

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
matUDIR = tempUINF./tempVELMAG;

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
rows = uint16(repmat(idxdve,1,m(1))) + uint16(tempm);

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
    airfoil = dlmread(strcat('airfoils/',cellAIRFOIL{pan},'.dat'),'', 1, 0);

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
% Apply direction to CDP
matCDPDIST = zeros(size(matUINF,1),3);
matCDPDIST(vecDVELE>0,:) = matUDIR(ledves,:).*vecCDPDIST;

end