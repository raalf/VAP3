function [matCENTER, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST, matNPVLST, matNTVLST, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM)

%   V0 - before fixing spanwise interp
%   V1 - fixed vertical panel (90deg dihedral)
%      - Reprogrammed the DVE interpolation method
%  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
%   V2 - Rework the Leading Edge Vector
%      - Caculate the Yaw angle of the DVE by assuming yaw=0 to rotate the xsi vector
%      - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
%      - Comptue LE Sweep
%      - Project TE to DVE, Rotate adn Comptue TE Sweep
%   V3 - Function overhaul for VAP2.0
%  3.5 - Modify non-planer VLST (Jan 6, 2017)
%        old matNPVLST will be not called matNTVLST to specify it holds dve
%        infomation of non-twisted wing
%      - new matNPVLST will now hold non-planer dve coordinates of
%        non-modified wing geometry specified in input file
%
% Fixed how DVEs matrix is converted from 2D grid to 1D array. 16/01/2016 (Alton)

% INPUT:
%   valPANELS - number of wing panels
%   matGEOM - 2 x 5 x valPANELS matrix, with (x,y,z) coords of edge points, and chord and twist at each edge
%   vecSYM - valPANELS x 1 vector of 0, 1, or 2 which denotes the panels with symmetry
%   vecN - valPANELS x 1 vector of spanwise elements per DVE
%   vecM - valPANELS x 1 vector of chordwise elements per DVE

% OUTPUT: (ALL OUTPUT ANGLES ARE IN RADIAN)
%   matCENTER - valNELE x 3 matrix of (x,y,z) locations of DVE control points
%   vecDVEHVSPN - valNELE x 1 vector of DVE half spans
%   vecDVEHVCRD - valNELE x 1 vector of DVE half chords
%   vecDVELESWP - valNELE x 1 vector of DVE leading edge sweep (radians)
%   vecDVEMCSWP - valNELE x 1 vector of DVE mid-chord sweep (radians)
%   vecDVETESWP - valNELE x 1 vector of DVE trailing-edge sweep (radians)
%   vecDVEROLL - valNELE x 1 vector of DVE roll angles (about x-axis) (radians)
%   vecDVEPITCH - valNELE x 1 vector of DVE pitch angles (about y-axis) (radians)
%   vecDVEYAW - valNELE x 1 vector of DVE yaw angles (about z-axis) (radians)
%   vecDVEAREA - valNELE x 1 vector of DVE area
%   vecDVENORM -  valNELE x 3 matrix of DVE normal vectors
%   matVLST - ? x 3 list of unique vertices, columns are (x,y,z) values
%   valNELE - total number of DVEs
%   matDVE - matrix of which DVE uses which vertices from the above list
%   matADJE - matADJE - ? x 3 adjacency matrix, where columns are: DVE | local edge | adjacent DVE
%   vecDVESYM - valNELE x 1 vector of which DVEs have symmetry on which edge (0 for no symmetry, 2 for local edge 2, 4 for local edge 4)
%   vecDVETIP - valNELE x 1 vector of which DVEs are at the wingtip. Similar format to vecDVESYM

% FUNCTIONS USED:
%   fcnPANELCORNERS
%   fcnPANEL2DVE
%   fcnGLOBSTAR

%% Preallocation
valNELE = sum(vecM.*vecN);


vecDVEPANEL   = nan(valNELE,1);

P1          = nan(valNELE,3);
P2          = nan(valNELE,3);
P3          = nan(valNELE,3);
P4          = nan(valNELE,3);

vecDVEWING  = nan(valNELE,1);




vecEnd      = cumsum(vecN.*vecM);




%% Assign Wing to Panel
panelEdges = reshape(permute(matGEOM,[1 3 2]),[],5);
[~,tempB,tempC] = unique(panelEdges,'rows','stable');
panelEdgesIdx = reshape(tempC,2,[])';
edge2wing = [(1:length(tempB))',nan(length(tempB),1)];
% Assign first edge to wing 1
edge2wing(1,2) = 1;

for n = 1:valPANELS
    curPanelEdges = panelEdgesIdx(n,:)';
    if max(edge2wing(curPanelEdges,2)) > 0
        wingIdx = max(edge2wing(curPanelEdges,2));
    else
        wingIdx = max(edge2wing(:,2))+1;
    end
    edge2wing(curPanelEdges,2)=wingIdx;
end
temp1 = reshape(edge2wing(tempC,2),2,[]);
panel2wing = temp1(1,:)';
clear tempB tempC temp1

%% Convert Panels to Corner Points to DVEs

for i = 1:valPANELS
    
    rchord = matGEOM(1,4,i); repsilon = deg2rad(matGEOM(1,5,i));
    tchord = matGEOM(2,4,i); tepsilon = deg2rad(matGEOM(2,5,i));
    rLE = matGEOM(1,1:3,i);
    tLE = matGEOM(2,1:3,i);
    
    % Read panel corners
    % For DVE generation. Twist angle is handled along with dihedral angle
    panel4corners = reshape(fcnPANELCORNERS(rLE,tLE,rchord,tchord,repsilon,tepsilon),3,4)';
    
    % fcnPANEL2DVE takes four corners of a panel and outputs vertices of non-planer DVEs
    [ CP, LE_Left, LE_Right, TE_Left, TE_Right ] = fcnPANEL2DVE( panel4corners, i, vecN, vecM );
    
    % Imaginary Wing for panel adjacencies. Twist of the panels are ignored
    % to ensure no gaps between panels on same wing.
    impanel4corners = reshape(fcnPANELCORNERS(rLE,tLE,rchord,tchord,0,0),3,4)';
    % fcnPANEL2DVE takes four corners of a panel and outputs vertices of non-planer DVEs
    [ ~, imLEL, imLER, imTEL, imTER ] = fcnPANEL2DVE( impanel4corners, i, vecN, vecM );
    
    % WRITE RESULTS
    count = vecN(i)*vecM(i);
    idxStart = vecEnd(i)-count+1;
    idxEnd = vecEnd(i);
    
    

    vecDVEPANEL(idxStart:idxEnd,:) = repmat(i,count,1);
    
    % Write DVE WING Index
    vecDVEWING(idxStart:idxEnd,:) = repmat(panel2wing(i),count,1);

    % Write DVE CENTER POINT Coordinates
    matCENTER(idxStart:idxEnd,:) = reshape(permute(CP, [2 1 3]),count,3);%reshape(CP(:),count,3);

    % Write non-planer DVE coordinates
    P1(idxStart:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count,3);
    P2(idxStart:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count,3);
    P3(idxStart:idxEnd,:) = reshape(permute(TE_Right, [2 1 3]),count,3);
    P4(idxStart:idxEnd,:) = reshape(permute(TE_Left, [2 1 3]),count,3);

    % Write Imeragary Wings
    imP1(idxStart:idxEnd,:) = reshape(permute(imLEL, [2 1 3]),count,3);
    imP2(idxStart:idxEnd,:) = reshape(permute(imLER, [2 1 3]),count,3);
    imP3(idxStart:idxEnd,:) = reshape(permute(imTER, [2 1 3]),count,3);
    imP4(idxStart:idxEnd,:) = reshape(permute(imTEL, [2 1 3]),count,3);
    
    clear LE_Left LE_Mid LE_Right TE_Right TE_Left ...
        imLEL imLER imTER imTEL ...
        idxStart idxEnd count
end


%% fcnDVECORNER2PARAM takes the corner and center points of each DVEs,
% computes the parameters and compiles the matVLST and matDVE

[ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, ...
    matVLST, matDVE, ~, idxVLST] = fcnDVECORNER2PARAM( matCENTER, P1, P2, P3, P4 );


%% Create nonplaner VLST
nonplanerVLST = [P1;P2;P3;P4];
matNPVLST = nonplanerVLST(idxVLST,:);



%% Solve ADJT DVE
% Grab the imaginary (no-twist) non-planer vertex list to avoid the gaps between DVEs
notwistnonplanerVLST = [imP1;imP2;imP3;imP4];
matNTVLST = notwistnonplanerVLST(idxVLST,:);


[ matADJE, vecDVESYM, vecDVETIP, vecDVELE, vecDVETE ] = fcnDVEADJT( imP1, imP2, imP3, imP4, valNELE, vecDVEPANEL, vecSYM );



end








