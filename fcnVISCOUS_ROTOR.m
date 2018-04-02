function [matDPDIST, vecDELNDIST] = fcnVISCOUS_ROTOR( valKINV, ...
    vecDVEHVCRD, vecN, vecM, vecDVELE, vecDVEPANEL, cellAIRFOIL, vecDISNORM, vecDVEAREA, matUINF, matVLST, matDVE, matWUINF)
% This function applies a viscous correction to rotors using look up tables
% and applies a high angle stall model.
%
% OUTPUT
%   matCDPDIST - CDP with direction accounted for
%   vecDELNDIST - The change in lift distribution due to stall

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

vecCDPDIST = nan(size(vecCNDIST)); % pre-allocate the array to store viscous drag results
vecCLMAXA   = nan(size(vecCNDIST));

% collect all the unique airfoils and load them
[uniqueAirfoil,~,idxAirfoil] = unique(cellAIRFOIL);
for k = 1:length(uniqueAirfoil)
    % Load airfoil .mat files
    try
        load(strcat('airfoils/',cellAIRFOIL{k},'.mat'));
        
    catch
        error('Error: Unable to locate airfoil file: %s.mat.', cellAIRFOIL{k});
    end
    
    Alpha  = reshape(pol(:,1,:),[],1);
    Cdp = reshape(pol(:,3,:),[],1);
    Re  = reshape(pol(:,8,:),[],1);
    
    
    %which rows of DVE belongs to the airfoil in this loop
    isCurrentAirfoil = idxAirfoil(idxpanel) == k;
    
    
    idxNans = isnan(Alpha) | isnan(Cdp) | isnan(Re);
    Alpha = Alpha(~idxNans);
    Cdp = Cdp(~idxNans);
    Re = Re(~idxNans);
    
    
    % Compare Re data range to panel Re
    if max(vecREDIST(isCurrentAirfoil)) > max(Re)
        disp('Re higher than airfoil Re data.')
    end
    
    if min(vecREDIST(isCurrentAirfoil)) < min(Re)
        disp('Re lower than airfoil Re data.')
    end
    
    
    % find CLmax for each row of dves
    [polarClmax, idxClmax] = max(pol(:,2,:));
    
    polarClmaxA = nan(size(polarClmax));
    for p = 1:length(polarClmax)
        polarClmaxA(p) = pol(idxClmax(:,:,p),1,p);
    end
    polarClmaxA = polarClmaxA(:);
    polarClmax = polarClmax(:);
    
    polarClmaxRe = unique(pol(:,8,:));
    polarClmaxRe = sort(polarClmaxRe(~isnan(polarClmaxRe)));
    %polarClmaxRe = reshape(pol(1,8,:),[],1);
    
    vecCLMAXA(isCurrentAirfoil) = interp1(polarClmaxRe,polarClmaxA,vecREDIST(isCurrentAirfoil),'linear');
    
    
    % Out of range Reynolds number index
    idxReOFR = (vecREDIST > max(Re) | vecREDIST < min(Re)) & isCurrentAirfoil;
    
    % Nearest extrap for out of range Reynolds number
    vecCLMAXA(idxReOFR) = interp1(polarClmaxRe,polarClmaxA,vecREDIST(idxReOFR),'nearest','extrap');
    
    % Check for stall and change the CL
    idxSTALL = (radtodeg(vecALPHAEFF) > vecCLMAXA) & isCurrentAirfoil;
    
    F = scatteredInterpolant(Re,Alpha,Cdp,'linear');
    
    if sum(idxSTALL) > 1
        disp('Airfoil sections have stalled.')      
        % Make apply stall model using empirical equations
        % cn = cd,90*(sin(alpha_eff))/(0.56+0.44sin(alpha_eff))
        % ct = cd,0*cos(alpha_eff)/2
        % cd = cn*sin(alpha_eff)+ct*cos(alpha_eff)
        % Note: cd,90 = 2
        % Find cd_0
        
        cd_0 = F(vecREDIST(idxSTALL),zeros(sum(idxSTALL),1));
        
        cn = 2.*sin(abs(vecALPHAEFF(idxSTALL)))./(0.56+0.44.*sin(abs(vecALPHAEFF(idxSTALL))));
        ct = cd_0.*cos(abs(vecALPHAEFF(idxSTALL)))./2;
        vecCNDIST0(idxSTALL)  = cn.*cos(abs(vecALPHAEFF(idxSTALL))) - ct.*sin(abs(vecALPHAEFF(idxSTALL)));
        
        vecCDPDIST(idxSTALL) = cn.*sin(abs(vecALPHAEFF(idxSTALL))) + ct.*cos(abs(vecALPHAEFF(idxSTALL)));
        
    end
    
    vecCDPDIST(idxSTALL==0) = F(vecREDIST(idxSTALL==0), radtodeg(vecALPHAEFF(idxSTALL==0)));
    
    clear pol foil
end

% Apply direction to CDP
matCDPDIST = zeros(size(matUINF,1),3);
matCDPDIST(vecDVELE>0,:) = matUDIR(ledves,:).*vecCDPDIST;

matDPDIST = 0.5*(sum(vecDVEAREA(rows),2).*(vecV.^2)).*matCDPDIST;
matDELCNDIST = (vecCNDIST0- vecCNDIST);
vecDELNDIST = 0.5*(sum(vecDVEAREA(rows),2).*(vecV.^2)).*matDELCNDIST;
end