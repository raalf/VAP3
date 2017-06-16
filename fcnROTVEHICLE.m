function [ matVLST, matCENTER, matFUSEGEOM] = fcnROTVEHICLE( matDVE, matVLST, matCENTER, valVEHICLES, vecDVEVEHICLE, matVEHORIG, matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS)
%FCNROTVEHICLE Summary of this function goes here
%   using for loop to rotate each vehicle with respect with their origin
%   position in global coordinates

for n = 1:valVEHICLES
    idxVLSTVEH = unique(matDVE(vecDVEVEHICLE==n,:));
    idxDVEVEH = vecDVEVEHICLE==n;
    valVLSTVEH = length(idxVLSTVEH);
    valDVEVEH = sum(idxDVEVEH);
    % glob to local translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) - repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) - repmat(matVEHORIG(n,:),valDVEVEH,1);
    
    % rotate
    dcm = angle2dcm(matVEHROT(n,1), matVEHROT(n,2), matVEHROT(n,3), 'XYZ');
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:)*dcm;
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:)*dcm;
    
    % local to global translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) + repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) + repmat(matVEHORIG(n,:),valDVEVEH,1);
end

for i = 1:length(vecFUSEVEHICLE)
    
    A = matFUSEAXIS(i,:);
    
    B = fcnSTARGLOB(matFUSEAXIS(i,:), matVEHROT(vecFUSEVEHICLE(i),1), matVEHROT(vecFUSEVEHICLE(i),2), matVEHROT(vecFUSEVEHICLE(i),3));
    
    if ~isequal(A,B)
        v = cross(A,B);
        ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + ssc + ssc^2*(1-dot(A,B))/(norm(v))^2;
    else
        R = [1 0 0; 0 1 0; 0 0 1];
    end
    
    X = matFUSEGEOM(:,:,1,i);
    Y = matFUSEGEOM(:,:,2,i);
    Z = matFUSEGEOM(:,:,3,i);
    
    temp=[X(:),Y(:),Z(:)]*R.';
    sz=size(X);
    Xrot=reshape(temp(:,1),sz);
    Yrot=reshape(temp(:,2),sz);
    Zrot=reshape(temp(:,3),sz);
    
    matFUSEGEOM(:,:,1,i) = Xrot + matVEHORIG(i,1);
    matFUSEGEOM(:,:,2,i) = Yrot + matVEHORIG(i,2);
    matFUSEGEOM(:,:,3,i) = Zrot + matVEHORIG(i,3);
    
    
end

