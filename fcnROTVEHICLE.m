function [matVLST, matCENTER,...
    matROTORHUBGLOB, matROTORAXIS, matNTVLST, vecBFRAME] = fcnROTVEHICLE( matDVE, matNPDVE, matVLST, ...
    matCENTER, valVEHICLES, vecDVEVEHICLE, matVEHORIG, ...
    matVEHROT, matROTORHUB, matROTORAXIS, vecROTORVEH, matNTVLST, vecBFRAME)
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
    idxNPVLSTVEH = unique(matNPDVE(vecDVEVEHICLE==n,:));
    idxDVEVEH = vecDVEVEHICLE==n;
    valVLSTVEH = length(idxVLSTVEH);
    valNPVLSTVEH = length(idxNPVLSTVEH);
    valDVEVEH = sum(idxDVEVEH);
    
    
    % glob to local translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) - repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matNTVLST(idxNPVLSTVEH,:) = matNTVLST(idxNPVLSTVEH,:) - repmat(matVEHORIG(n,:),valNPVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) - repmat(matVEHORIG(n,:),valDVEVEH,1);
    
    % rotate
    %     dcm = angle2dcm(matVEHROT(n,1), matVEHROT(n,2), matVEHROT(n,3), 'XYZ');
    dcm = angle2dcm(matVEHROT(n,3), matVEHROT(n,1), matVEHROT(n,2), 'ZXY');
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:)*dcm;
    matNTVLST(idxNPVLSTVEH,:) = matNTVLST(idxNPVLSTVEH,:)*dcm;
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:)*dcm;
    
    dcm = angle2dcm(matVEHROT(n,3), matVEHROT(n,1), -matVEHROT(n,2), 'ZXY');
    vecBFRAME = dcm*vecBFRAME;
    
    
    if ~isempty(matROTORHUB)
        % rotate vecROTORVEH for later reference
        matROTORHUBGLOB(vecROTORVEH==n,:) = matROTORHUB(vecROTORVEH==n,:)*dcm;
    end
    
    % local to global translation
    matVLST(idxVLSTVEH,:) = matVLST(idxVLSTVEH,:) + repmat(matVEHORIG(n,:),valVLSTVEH,1);
    matNTVLST(idxNPVLSTVEH,:) = matNTVLST(idxNPVLSTVEH,:) + repmat(matVEHORIG(n,:),valNPVLSTVEH,1);
    matCENTER(idxDVEVEH,:) = matCENTER(idxDVEVEH,:) + repmat(matVEHORIG(n,:),valDVEVEH,1);
end

