function [valCL, valCD, valCM, valPREQ, valLD, valVINF, vecCMDIST, temp_CN, vecCDPDIST, temp_dist] = fcnVISCOUS_WING(valCL, valCDI, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
    matCENTER, vecDVEHVCRD, cellAIRFOIL, flagPRINT, vecSYM, ...
    valINTERF, vecDVEROLL, matUINF, matWUINF, matDVE, matVLST, valVEHVINF, fixed_lift, valVEHWEIGHT, vecCMAC, valFUSEFPA)

warning off

% % % Calculate chordline direction at midspan of each dve
% % avgle = (matVLST(matDVE(:,1),:)+matVLST(matDVE(:,2),:))./2;
% % avgte = (matVLST(matDVE(:,3),:)+matVLST(matDVE(:,4),:))./2;
% % tempDIF = avgte - avgle;
% % matCRDLINE = (tempDIF)./(repmat(sqrt(sum(tempDIF.^2,2)),[1,3]));
% Calculate velocity with induced effects
%vecV = dot(matUINF(vecDVEWING>0) + matWUINF(vecDVEWING>0), matCRDLINE(vecDVEWING>0),2);
valVINF = nan;
if fixed_lift ~= 1
    temp  = sqrt(matUINF(:,1).^2 + matUINF(:,2).^2 + matUINF(:,3).^2);
    temp1 = dot(matWUINF, matUINF./temp, 2);
    vecV = temp1 + temp;
    
    % Compute dynamic pressure
    q_infandind = ((vecV.^2)*valDENSITY)/2; % With both freestream and induced velocities
    q_inf = ((valVEHVINF^2)*valDENSITY)/2; % With only freestream velocities
    valVINF = mean(vecV,1);
else
    q_inf = valVEHWEIGHT./(valCL.*valAREA);
    q_infandind = repmat(q_inf, size(matUINF,1), 1);
    valVINF = sqrt(2.*q_inf./valDENSITY);
    
    if flagPRINT == 1
        disp(['Using janky fixed-lift analysis - VINF = ', num2str(valVINF)]);
    end
    
    temp  = sqrt(matUINF(:,1).^2 + matUINF(:,2).^2 + matUINF(:,3).^2);
    temp1 = dot(matWUINF, matUINF./temp, 2);
    vecV = temp1 + temp;
end

% Calculate induced drag as a force
di = valCDI*valAREA*q_inf;

% Summing freestream and induced forces of each DVE
vecDVECN = (vecDVENFREE + vecDVENIND);

% [ledves, ~, ~] = find(vecDVELE > 0);
[ledves, ~, ~] = find(vecDVELE > 0 & vecDVEWING > 0);
lepanels = vecDVEPANEL(ledves);

vecCNDIST    = nan(size(ledves,1),1);
matXYZDIST   = nan(size(ledves,1),3);
vecLEDVEDIST = nan(size(ledves,1),1);
vecREDIST    = nan(size(ledves,1),1);
vecAREADIST  = nan(size(ledves,1),1);
vecCDPDIST   = nan(size(ledves,1),1); % pre-allocate the array to store viscous drag
vecCMDIST    = nan(size(ledves,1),1);
vecCLMAX     = nan(size(ledves,1),1);
dprofPerWing = nan(max(vecDVEWING),1);
LPerWing     = nan(max(vecDVEWING),1);

for i = 1:max(vecDVEWING)
    
    %% Getting the CL, CY, CN distribution
    isCurWing = vecDVEWING(ledves) == i;
    
    idxdve = uint16(ledves(isCurWing));
    idxpanel = lepanels(isCurWing);
    
    
    m = vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    
    m = m(1);
    
    % Matrix of how much we need to add to an index to get the next chordwise element
    % It is done this way because n can be different for each panel. Unlike in the wake,
    % we can't just add a constant value to get to the same spanwise location in the next
    % row of elements
    tempm = repmat(vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);
    
    rows = repmat(idxdve,1,m) + uint16(tempm);
    
    % Note this CN is non-dimensionalized with Vinf + Vind
    vecCNDIST(isCurWing) = (sum(vecDVECN(rows),2).*2)./((mean(vecV(rows),2).^2).*sum(vecDVEAREA(rows),2));
    
    % The average coordinates for this row of elements
    matXYZDIST(isCurWing,:) = mean(permute(reshape(matCENTER(rows,:)',3,[],m),[2 1 3]),3);
    
    % The leading edge DVE for the distribution
    vecLEDVEDIST(isCurWing) = idxdve;
    
    %% Wing/horizontal stabilizer lift and drag
    % Note that Re is compute with Vinf + Vind
    vecREDIST(isCurWing)   = mean(vecV(rows),2).*2.*sum(vecDVEHVCRD(rows),2)./valKINV;
    % If fixed-lift is enabled, UINF is defined as unity.
    % Re values have to be scaled to have correct viscous results.
    if fixed_lift == 1
        vecREDIST(isCurWing) = vecREDIST(isCurWing)*valVINF;
    end
    
    vecAREADIST(isCurWing) = sum(vecDVEAREA(rows),2);
    
    vecCNDIST0 = vecCNDIST;
    
    % collect all the unique airfoils and load them
    [uniqueAirfoil,~,idxAirfoil] = unique(cellAIRFOIL, 'stable');
    for k = 1:length(uniqueAirfoil)
        % Load airfoil .mat files
        try
            % only load variable 'pol' to avoid variable conflict.
            load(strcat('airfoils/',uniqueAirfoil{k},'.mat'),'pol');
            
        catch
            error('Error: Unable to locate airfoil file: %s.mat.', uniqueAirfoil{k});
        end
        
        Cl  = reshape(pol(:,2,:),[],1);
        Cdp = reshape(pol(:,3,:),[],1);
        Cm = reshape(pol(:,5,:),[],1);
        Re  = reshape(pol(:,8,:),[],1);
        
        %which rows of DVE belongs to the airfoil in this loop
        isCurrentAirfoil = isCurWing & idxAirfoil(lepanels) == k;
        
        idxNans = isnan(Cl) | isnan(Cdp) | isnan(Re);
        Cl = Cl(~idxNans);
        Cdp = Cdp(~idxNans);
        Cm = Cm(~idxNans);
        Re = Re(~idxNans);
        
        % Compare Re data range to panel Re
        %         if max(vecREDIST(isCurrentAirfoil)) > max(Re)
        %             disp('fcnVISCOUS_WING: Re higher than airfoil Re data.')
        %         end
        %
        %         if min(vecREDIST(isCurrentAirfoil)) < min(Re)
        %             disp('fcnVISCOUS_WING: Re lower than airfoil Re data.')
        %         end
        
        % find CLmax for each row of dves
        polarClmax = max(pol(:,2,:));
        polarClmax = polarClmax(:);
        
        polarClmaxRe = unique(pol(:,8,:));
        polarClmaxRe = sort(polarClmaxRe(~isnan(polarClmaxRe)));
        %polarClmaxRe = reshape(pol(1,8,:),[],1);
        
        % Remove NANs in polarClmax and index polarClmaxRe
        polarClmax   = polarClmax(~isnan(polarClmax));
        polarClmaxRe = polarClmaxRe(~isnan(polarClmax));
        
        vecCLMAX(isCurrentAirfoil) = interp1(polarClmaxRe,polarClmax,vecREDIST(isCurrentAirfoil),'linear');
        
        if any(isCurrentAirfoil)
            idxNans = isnan(Cl) | isnan(Cdp) | isnan(Re);
            Cl = Cl(~idxNans);
            Cdp = Cdp(~idxNans);
            Re = Re(~idxNans);
            
            % Compare Re data range to panel Re
            if (max(vecREDIST(isCurrentAirfoil)) > max(Re)) && flagPRINT == 1
%                 fprintf('fcnVISCOUS_WING: Re higher than airfoil Re data. Re = %i, (%s)\n', max(vecREDIST(isCurrentAirfoil)), uniqueAirfoil{k})
            end
            
            if (min(vecREDIST(isCurrentAirfoil)) < min(Re)) && flagPRINT == 1
%                 fprintf('fcnVISCOUS_WING: Re lower than airfoil Re data. Re = %i, (%s)\n',  min(vecREDIST(isCurrentAirfoil)), uniqueAirfoil{k})
            end
            
            % find CLmax for each row of dves
            polarClmax = max(pol(:,2,:));
            polarClmax = polarClmax(:);
            
            polarClmaxRe = unique(pol(:,8,:));
            polarClmaxRe = sort(polarClmaxRe(~isnan(polarClmaxRe)));
            %polarClmaxRe = reshape(pol(1,8,:),[],1);
            
            % Remove NANs in polarClmax and index polarClmaxRe
            polarClmax   = polarClmax(~isnan(polarClmax));
            polarClmaxRe = polarClmaxRe(~isnan(polarClmax));
            
            vecCLMAX(isCurrentAirfoil) = interp1(polarClmaxRe,polarClmax,vecREDIST(isCurrentAirfoil),'linear');
            
            % Out of range Reynolds number index
            idxReOFR = (vecREDIST > max(Re) | vecREDIST < min(Re)) & isCurrentAirfoil;
            % Nearest extrap for out of range Reynolds number
            vecCLMAX(idxReOFR) = interp1(polarClmaxRe,polarClmax,vecREDIST(idxReOFR),'nearest','extrap');
            
            % Check for stall and change the CL
            idxSTALL = (vecCNDIST > vecCLMAX) & isCurrentAirfoil;
            vecCNDIST0(idxSTALL) = vecCLMAX(idxSTALL)*0.825;
            if sum(idxSTALL) > 1 && flagPRINT == 1
                fprintf('fcnVISCOUS_WING: Airfoil sections have stalled. (%s)\n', uniqueAirfoil{k})
            end
            
            F = scatteredInterpolant(Re,Cl,Cdp,'linear','nearest');
            vecCDPDIST(isCurrentAirfoil) = F(vecREDIST(isCurrentAirfoil), vecCNDIST(isCurrentAirfoil));
            clear pol foil
        end
        
        F = scatteredInterpolant(Re,Cl,Cdp,'linear','nearest');
        vecCDPDIST(isCurrentAirfoil) = F(vecREDIST(isCurrentAirfoil), vecCNDIST(isCurrentAirfoil));
        
        F = scatteredInterpolant(Re,Cl,Cm,'linear','nearest');
        vecCMDIST(isCurrentAirfoil) = F(vecREDIST(isCurrentAirfoil), vecCNDIST(isCurrentAirfoil));
        clear pol foil
    end
    % CN in terms of Vinf instead of Vinf + Vind
    vecCNDIST(isCurWing) = vecCNDIST0(isCurWing).*(mean(vecV(rows),2).^2)/(valVEHVINF^2);
    
%     temp_dist(i).LDIST = vecCNDIST(isCurWing).*cos(vecDVEROLL(vecLEDVEDIST(isCurWing))).*q_inf.*0.395%vecAREADIST(isCurWing);
    temp_dist(i).LDIST = vecCNDIST(isCurWing).*cos(vecDVEROLL(vecLEDVEDIST(isCurWing))).*q_inf.*vecAREADIST(isCurWing);
    temp_dist(i).DPDIST = vecCDPDIST(isCurWing).*mean(q_infandind(rows),2).*vecAREADIST(isCurWing);
    
    % Lift force per wing
    LPerWing(i) = sum(vecCNDIST(isCurWing).*cos(vecDVEROLL(vecLEDVEDIST(isCurWing))).*q_inf.*vecAREADIST(isCurWing));
    
    % dimensionalize in terms of both Vinf and Vind
    dprofPerWing(i) = sum(vecCDPDIST(isCurWing).*mean(q_infandind(rows),2).*vecAREADIST(isCurWing));
    pitchmomPerWing(i) = sum(vecCMDIST(isCurWing).*mean(q_infandind(rows),2).*vecAREADIST(isCurWing)*vecCMAC);
    temp_CN = vecCNDIST;

end



% create temporary variables to determine which panel belongs to which wing
% in order to work out the symmetry per wing

temp1 = unique([vecDVEPANEL,vecDVEWING],'rows');
vecPanelWingNumber(temp1(:,1),1) = temp1(:,2);
vecSYMWING = false(max(vecPanelWingNumber),1);
for j = 1:max(vecPanelWingNumber)
    vecSYMWING(j,1) = any(vecSYM(vecPanelWingNumber == j));
end

% Multiply vecSYMIWING profile drag by 2
dprofPerWing(vecSYMWING) = dprofPerWing(vecSYMWING)*2;
pitchmomPerWing(vecSYMWING) = pitchmomPerWing(vecSYMWING)*2;
LPerWing(vecSYMWING) = LPerWing(vecSYMWING)*2;

valCL = sum(LPerWing)/(q_inf*valAREA);

% sum profile drag per wing
dprof = sum(dprofPerWing);

%% Total Drag

dtot = di + dprof + valFUSEFPA*q_inf;

dint = dtot*(valINTERF/100);

dtot = dtot + dint;

valCD = dtot/(q_inf*valAREA);

valCM = sum(pitchmomPerWing)/(q_inf*valAREA*vecCMAC);

% %% Adjusting CL for stall
% valCL = sum(vecCNDIST.*vecAREADIST.*cos(vecDVEROLL(vecLEDVEDIST)))/valAREA*2;

%% Final calculations
valLD = valCL./valCD;
if fixed_lift == 1
    valPREQ = dtot.*valVINF;
else
    valPREQ = dtot.*valVEHVINF;
end

end

