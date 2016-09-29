function [matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER, matVLST, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecSYM)

% Flow tangency is to be enforced at all control points on the surface HDVEs
% In the D-Matrix, dot (a,b,c) of our influencing HDVE with the normal of the point we are influencing on

%% Adding king kong conditions to bottom 1/3 of D-matrix

% Points we are influencing:
fpg = matCENTER;

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);

fpg = repmat(fpg,valNELE,1);

% DVE type 0 is a surface element
dvetype = zeros(length(dvenum),1);

[a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, vecK, matDVE, matVLST, vecDVEHVSPN, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecSYM);

% List of normals we are to dot the above with
normals = repmat(matDVENORM,valNELE,1); % Repeated so we can dot all at once

% Dotting a, b, c with the normals of the field points
temp60 = [dot(a,normals,2) dot(b,normals,2) dot(c,normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*3);
king_kong(rows,:) = reshape(permute(reshape(temp60',3,[],valNELE),[2 1 3]),[],3*valNELE,1);

matD = [matD; king_kong];

end

