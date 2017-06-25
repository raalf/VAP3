function [ matUINF ] = fcnINITUINF( matCENTER, matVEHUVW, matVEHROT, vecDVEVEHICLE, ...
    vecDVEROTOR, vecROTORVEH, matVEHORIG, matROTORHUBGLOB, matROTORAXIS, vecROTORRPM )

%FCNINITUINF Summary of this function goes here
%   initializes DVE UINF by adding vehicle UVW and calculating rotor UVW
matUINF = -matVEHUVW(vecDVEVEHICLE,:);
valROTORS = length(vecROTORVEH);
vecROTORRADPS = vecROTORRPM.*2.*pi./60;

for n = 1:valROTORS
    idxDVEROTOR = vecDVEROTOR==n;
    tempROTORCENTER = matCENTER(idxDVEROTOR,:);
    tempROTORCENTER = tempROTORCENTER - matROTORHUBGLOB(n,:) - matVEHORIG(vecROTORVEH(n),:);
    
    % transform rotor from global to hub plane
    tempROTORCENTER = tempROTORCENTER / angle2dcm(matVEHROT(vecROTORVEH(n),3),matVEHROT(vecROTORVEH(n),1),matVEHROT(vecROTORVEH(n),2),'ZXY');    
    
    % transform rotor from hub plane to xy plane
    tempROTORCENTER = tempROTORCENTER * quat2dcm(axang2quat(vrrotvec([0 0 1],matROTORAXIS(n,:))));

    % timestep rotor in local XY hub plane
    tempROTORUINF = cross(repmat([0,0,-vecROTORRADPS(n)],length(tempROTORCENTER(:,1)),1),tempROTORCENTER);    
    
    % transform rotor from xy plane to hub plane
	tempROTORUINF = tempROTORUINF * quat2dcm(axang2quat(vrrotvec(matROTORAXIS(n,:),[0 0 1])));
    
    % transform rotor from hub plane to global
    tempROTORUINF = tempROTORUINF * angle2dcm(matVEHROT(vecROTORVEH(n),3),matVEHROT(vecROTORVEH(n),1),matVEHROT(vecROTORVEH(n),2),'ZXY');
    
    % write rotated rotor to matVLST
    matUINF(idxDVEROTOR,:) = tempROTORUINF;
end


end

