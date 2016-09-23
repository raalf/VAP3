function [vecDVECTLPT,vecDVEHVSPN,vecDVEHVCRD,vecDVELESWP,vecDVEMCSWP,vecDVETESWP,vecDVEROLL,vecDVEPITCH,vecDVEYAW,vecDVEAREA,vecDVENORM,matVLST,matVIDX] = fcnGenerateDVEs(valPANELS, matGEOM, vecN, vecM)
% % %   V0 - before fixing spanwise interp
% % %   V1 - fixed vertical panel (90deg dihedral)
% % %      - Reprogrammed the DVE interpolation method
% % %  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
% % %   V2 - Rework the Leading Edge Vector
% % %      - Caculate the Yaw angle of the DVE by assuming yaw=0 to rotate the xsi vector
% % %      - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
% % %      - Comptue LE Sweep
% % %      - Project TE to DVE, Rotate adn Comptue TE Sweep
% % %   V3 - Function overhaul for VAP2.0
% %
% % % Fixed how DVEs matrix is converted from 2D grid to 1D array. 16/01/2016 (Alton)

% INPUT:
% valPANELS
% matGEOM
% vecN
% vecM
 
% OUTPUT:
% vecDVECTLPT
% vecDVEHVSPN
% vecDVEHVCRD
% vecDVELESWP
% vecDVEMCSWP
% vecDVETESWP
% vecDVEROLL
% vecDVEPITCH
% vecDVEYAW
% vecDVEAREA
% vecDVENORM
% matVLST
% matVIDX

dveCount = sum(vecM.*vecN);
vecDVECTLPT = nan(dveCount,3);
vecDVEHVSPN = nan(dveCount,1);
vecDVEHVCRD = nan(dveCount,1);
vecDVELESWP = nan(dveCount,1);
vecDVEMCSWP = nan(dveCount,1);
vecDVETESWP = nan(dveCount,1);
vecDVEROLL  = nan(dveCount,1);
vecDVEPITCH = nan(dveCount,1);
vecDVEYAW   = nan(dveCount,1);
vecDVEAREA  = nan(dveCount,1);
vecDVENORM  = nan(dveCount,3);

LECoordL    = nan(dveCount,3);
LECoordR    = nan(dveCount,3);
TECoordR    = nan(dveCount,3);
TECoordL    = nan(dveCount,3);

vecEnd = cumsum(vecN.*vecM);


for i = 1:valPANELS;

    rchord = matGEOM(1,4,i); repsilon = deg2rad(matGEOM(1,5,i));
    tchord = matGEOM(2,4,i); tepsilon = deg2rad(matGEOM(2,5,i));
    rLE = matGEOM(1,1:3,i);
    tLE = matGEOM(2,1:3,i);

    % Read panel corners
    panel = reshape(fcnPanelCorners(rLE,tLE,rchord,tchord,repsilon,tepsilon),3,4)';
    panelX = [panel([1;4],1),panel([2;3],1)];
    panelY = [panel([1;4],2),panel([2;3],2)];
    panelZ = [panel([1;4],3),panel([2;3],3)];
    %     panelTE = [panel(end,:); panel(end-1,:)];
    
    % NChordwise and NSpanwise
    % Generate extra points to find control point
    N = vecN(i)*2;
    M = vecM(i)*2;
    
    % split Root Chord to chordwise elements
    chordX(:,1) = linspace(panelX(1,1),panelX(2,1),M+1)';
    chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
    chordZ(:,1) = linspace(panelZ(1,1),panelZ(2,1),M+1)';
    % split Tip Chord to chordwise elements
    chordX(:,2) = linspace(panelX(1,2),panelX(2,2),M+1)';
    chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
    chordZ(:,2) = linspace(panelZ(1,2),panelZ(2,2),M+1)';
    
    % linspace spanwise elemenets at each chordwise station
    
    % Preallocate the memories for point matrices
    PX2 = zeros(M+1,N+1);
    PY2 = zeros(M+1,N+1);
    PZ2 = zeros(M+1,N+1);
    
    for j = 1:M+1
        if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
            spanwise = linspace(chordY(1,1),chordY(1,2),N+1);
            chordbase = chordY(1,:);
        else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
            spanwise = linspace(chordZ(1,1),chordZ(1,2),N+1);
            chordbase = chordZ(1,:);
        end
        PX2(j,:) = interp1(chordbase,chordX(j,:),spanwise); % X coordinates of all points
        PY2(j,:) = interp1(chordbase,chordY(j,:),spanwise); % Y coordinates of all points
        PZ2(j,:) = interp1(chordbase,chordZ(j,:),spanwise); % Z coordinates of all points
    end
    
    % DVE Parameters Calculation
    % Calculate Control Points, stored in 3D matrix
    CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],vecM(i),vecN(i),3);
    %     CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],vecM(i),vecN(i),3);
    LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],vecM(i),vecN(i),3);
    TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],vecM(i),vecN(i),3);
    TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],vecM(i),vecN(i),3);
    
    % Filter Leading Edge Point
    LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],vecM(i),vecN(i),3);
    LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],vecM(i),vecN(i),3);
    
    % Create eta vector for full leading edge
    % Non-normalized
    LE_vec = LE_Right - LE_Left;
    
    % Create half chord xsi vector
    % Non-normalized
    xsi_vec = LE_Mid - CP;
    
    
    %
    DVE_norm = fcnNORM3D(cross(LE_vec, xsi_vec, 3));
    
    
    
    % Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
    % atan2d is used here
    % roll(nu) right wing up positive
    nu = -atan2d(DVE_norm(:,:,2),DVE_norm(:,:,3));
    
    % Pitch in Degrees
    % arcsin ( X component of DVE normal vector )
    epsilon = asind(DVE_norm(:,:,1));
    
    % Yaw in Degrees
    % xsi in local with roll picth, yaw set to zero.. but WHY?
    xsi_local = fcnGLOBSTAR3D( xsi_vec,nu,epsilon,zeros(vecM(i),vecN(i)) );
    % Magnitude of half chord vector
    xsi = (xsi_local(:,:,1).^2+xsi_local(:,:,2).^2+xsi_local(:,:,3).^2).^0.5;
    psi = atand(xsi_local(:,:,2)./xsi_local(:,:,1));
    
    % Find eta. bring non-normalized LE_vec to local and half the Y component
    LE_vec_local = fcnGLOBSTAR3D( LE_vec,nu,epsilon,psi);
    eta = LE_vec_local(:,:,2)./2;
    
    % Find Leading Edge Sweep
    % arctan(LE X local component/ LE Y local component)
    phi_LE = atand(LE_vec_local(:,:,1)./LE_vec_local(:,:,2));
    
    
    
    % Find Trailing Edge Sweep
    % Project TE Points onto DVE plane
    % (TE_Left / TE_Right) (CP)                   (DVE_norm)
    % q(x,y,z) TE point | p(a,b,c) Control Point | n(d,e,f) DVE normal
    % q_proj = q - dot(q-p,n)*n
    TE_Left_proj = TE_Left-repmat(dot(TE_Left-CP,DVE_norm,3),1,1,3).*DVE_norm;
    TE_Right_proj = TE_Right-repmat(dot(TE_Right-CP,DVE_norm,3),1,1,3).*DVE_norm;
    TE_vec_proj = TE_Right_proj - TE_Left_proj;
    
    % Rotate the Projected TE on DVE to local reference frame
    % arctan(Projected TE local X component/Projected TE local Y component)
    TE_vec_proj_local = fcnGLOBSTAR3D( TE_vec_proj,nu,epsilon,psi );
    phi_TE = atand(TE_vec_proj_local(:,:,1)./TE_vec_proj_local(:,:,2));
    
    
    % Compute DVE Mid-chord Sweep
    % Average of LE and TE Sweep
    phi_MID = (phi_LE+phi_TE)./2;
    
    % Calculating Area
    Area = eta.*xsi.*4;
    
    
    
    count = vecN(i)*vecM(i);
    idxStart = vecEnd(i)-count+1;
    idxEnd = vecEnd(i);
    
    
    vecDVEAREA(idxStart:idxEnd) = reshape(Area',count,1);%Area(:);
    vecDVEHVSPN(idxStart:idxEnd) = reshape(eta',count,1);%eta(:);
    vecDVEHVCRD(idxStart:idxEnd) = reshape(xsi',count,1);%xsi(:);
    vecDVEROLL(idxStart:idxEnd) = reshape(nu',count,1);%nu(:);
    vecDVEPITCH(idxStart:idxEnd) = reshape(epsilon',count,1);%epsilon(:);
    vecDVEYAW(idxStart:idxEnd) = reshape(psi',count,1);%psi(:);
    vecDVECTLPT(idxStart:idxEnd,:) = reshape(permute(CP, [2 1 3]),count,3);%reshape(CP(:),count,3);
    vecDVELESWP(idxStart:idxEnd) = reshape(phi_LE',count,1);%phi_LE(:);
    vecDVEMCSWP(idxStart:idxEnd) = reshape(phi_MID',count,1);%phi_MID(:);
    vecDVETESWP(idxStart:idxEnd) = reshape(phi_TE',count,1);%phi_TE(:);
    vecDVENORM(idxStart:idxEnd,:) = reshape(permute(DVE_norm, [2 1 3]),count,3);%reshape(DVE_norm(:),count,3);
    
    
    LECoordL(idxStart:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
    LECoordR(idxStart:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);
    TECoordR(idxStart:idxEnd,:) = reshape(permute(TE_Right_proj, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);
    TECoordL(idxStart:idxEnd,:) = reshape(permute(TE_Left_proj, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
    
    clear chordX chordY chordZ panelX panelY panelZ m n debug ...
        X1 X2 Y1 Y2 Z1 Z2 AX1 AZ1 AX2 AZ2 C1 C2 e2 ...
        chordbase count PX2 PY2 PZ2 rLE ...
        panel paneldata i j M N ...
        halfchord halfspan
    
end

verticeList = [LECoordL;LECoordR;TECoordR;TECoordL];
[matVLST,~,matVIDX] = unique(verticeList,'rows');
matVIDX = reshape(matVIDX,dveCount,4);



end








