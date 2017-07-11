function [ matVLST, matCENTER, matFUSEGEOM, ...
    matROTORHUBGLOB, matROTORAXIS] = fcnROTVEHICLE( matDVE, matVLST, ...
    matCENTER, valVEHICLES, vecDVEVEHICLE, matVEHORIG, ...
    matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS, ...
    matROTORHUB, matROTORAXIS, vecROTORVEH)
%FCNROTVEHICLE Summary of this function goes here
%   using for loop to rotate each vehicle with respect with their origin
%   position in global coordinates

if ~isempty(matROTORHUB)
    matROTORHUBGLOB = nan(length(matROTORHUB(:,1)),3);
else
    matROTORHUBGLOB = [];
end

for n = 1:valVEHICLES
    idxVLSTVEH = unique(matDVE(vecDVEVEHICLE==n,:));
    idxDVEVEH = vecDVEVEHICLE==n;
    valVLSTVEH = length(idxVLSTVEH);
    valDVEVEH = sum(idxDVEVEH);
    
    
    % glob to local translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) - repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) - repmat(matVEHORIG(n,:),valDVEVEH,1);
    
    % rotate
%     dcm = angle2dcm(matVEHROT(n,1), matVEHROT(n,2), matVEHROT(n,3), 'XYZ');
    dcm = angle2dcm(matVEHROT(n,3), matVEHROT(n,1), matVEHROT(n,2), 'ZXY');
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:)*dcm;
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:)*dcm;
    
    if ~isempty(matROTORHUB)
        % rotate vecROTORVEH for later reference
        matROTORHUBGLOB(vecROTORVEH==n,:) = matROTORHUB(vecROTORVEH==n,:)*dcm;
    end
    
    % local to global translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) + repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) + repmat(matVEHORIG(n,:),valDVEVEH,1); 
    
    for i = 1:length(vecFUSEVEHICLE)
        if vecFUSEVEHICLE(i) == n
                X = matFUSEGEOM(:,:,1,i);
                Y = matFUSEGEOM(:,:,2,i);
                Z = matFUSEGEOM(:,:,3,i);
                temp=[X(:),Y(:),Z(:)]*dcm;
                sz=size(X);
                Xrot=reshape(temp(:,1),sz);
                Yrot=reshape(temp(:,2),sz);
                Zrot=reshape(temp(:,3),sz);
                matFUSEGEOM(:,:,1,i) = Xrot + matVEHORIG(i,1);
                matFUSEGEOM(:,:,2,i) = Yrot + matVEHORIG(i,2);
                matFUSEGEOM(:,:,3,i) = Zrot + matVEHORIG(i,3);
        end
    end
end

