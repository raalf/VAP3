function [matFUSEGEOM] = fcnCREATEFUSE(matSECTIONFUSELAGE, vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE)
%FCNCREATEFUSE Summary of this function goes here
%   Detailed explanation goes here
% hFig20 = figure(20);
% clf(20);

matFUSEGEOM = zeros(51, 21, 3, length(vecFUSESECTIONS));
for i = 1:length(vecFUSESECTIONS)
    
    [xData, yData] = prepareCurveData(matFGEOM(matSECTIONFUSELAGE == i,1), matFGEOM(matSECTIONFUSELAGE == i,2)./2);
    
    ft = fittype( 'poly3' );
    [fitresult, gof] = fit( xData, yData, ft );
    
    fuse_start = min(matFGEOM(matSECTIONFUSELAGE==1,1));
    fuse_end = max(matFGEOM(matSECTIONFUSELAGE==1,1));
    fuse_len = fuse_end - fuse_start;
    
    t = fuse_start:fuse_len/50:fuse_end;
    
    [X,Y,Z] = cylinder(fitresult(t));
    Z = Z.*fuse_len;
    
    A = [0 0 1];
    B = matFUSEAXIS(i,:);
    
    %     B = fcnSTARGLOB(matFUSEAXIS(i,:), matVEHROT(vecFUSEVEHICLE(i),1), matVEHROT(vecFUSEVEHICLE(i),2), matVEHROT(vecFUSEVEHICLE(i),3));
    
    v = cross(A,B);
    ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + ssc + ssc^2*(1-dot(A,B))/(norm(v))^2;
    
    temp=[X(:),Y(:),Z(:)]*R.';
    sz=size(X);
    Xrot=reshape(temp(:,1),sz);
    Yrot=reshape(temp(:,2),sz);
    Zrot=reshape(temp(:,3),sz);
    
    matFUSEGEOM(:,:,1,i) = Xrot + matFUSEORIG(i,1);
    matFUSEGEOM(:,:,2,i) = Yrot + matFUSEORIG(i,2);
    matFUSEGEOM(:,:,3,i) = Zrot + matFUSEORIG(i,3);
    
    
    %     surf(matFUSEGEOM(:,:,1,i), matFUSEGEOM(:,:,2,i), matFUSEGEOM(:,:,3,i))
    %     hold on
    
end
% hold off
% axis equal
% xlabel('X-Dir','FontSize',15);
% ylabel('Y-Dir','FontSize',15);
% zlabel('Z-Dir','FontSize',15);

end

