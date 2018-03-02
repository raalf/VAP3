function [ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, matVLST, matDVE, valNELE, idxVLST] = ...
    fcnDVECORNER2PARAM( matCENTER, P1, P2, P3, P4, vecDVEWING )

%FCNDVEPOINT2PARAM takes the corner and center points of each DVEs,computes the parameters and compiles the matVLST and matDVE
%   Detailed explanation goes here



valNELE = length(matCENTER(:,1)); 


% Create eta vector for full leading edge
% Non-normalized
% (old method) LE_vec = LE_Right - LE_Left;
leVec = P2-P1;

% Create half chord xsi vector
% Non-normalized
% (old method) xsi_vec = LE_Mid - CP;
P12 = (P1+P2)./2;
xsiVec = P12-matCENTER;



tempM = cross(leVec, xsiVec, 2);
tempM_length = repmat(((tempM(:,1).^2+tempM(:,2).^2+tempM(:,3).^2).^0.5),1,3);
matDVENORM = tempM./tempM_length;



% Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
% atan2d is used here
% roll(nu) right wing up positive
% (old method) nu = -atan2(DVE_norm(:,:,2),DVE_norm(:,:,3));
vecDVEROLL = -atan2(matDVENORM(:,2),matDVENORM(:,3));


% Pitch in Degrees
% arcsin ( X component of DVE normal vector )
% (old method) epsilon = asin(DVE_norm(:,:,1));
vecDVEPITCH = asin(matDVENORM(:,1));


% Yaw in Degrees
% xsi in local with roll picth, yaw set to zero.. but WHY?
% (old method) xsi_local = fcnGLOBSTAR3D( xsi_vec,nu,epsilon,zeros(vecM(i),vecN(i)) );
xsiLocal = fcnGLOBSTAR(xsiVec,vecDVEROLL,vecDVEPITCH,zeros(valNELE,1));
% % Magnitude of half chord vector
% (old method) xsi = (xsi_local(:,:,1).^2+xsi_local(:,:,2).^2+xsi_local(:,:,3).^2).^0.5;
vecDVEHVCRD = (xsiLocal(:,1).^2+xsiLocal(:,2).^2+xsiLocal(:,3).^2).^0.5;
% (old method) psi = atan(xsi_local(:,:,2)./xsi_local(:,:,1));
vecDVEYAW = atan2(-xsiLocal(:,2), -xsiLocal(:,1));
% 
% % Find eta. bring non-normalized LE_vec to local and half the Y component
% (old method) LE_vec_local = fcnGLOBSTAR3D( LE_vec,nu,epsilon,psi);
leVecLocal = fcnGLOBSTAR(leVec,vecDVEROLL,vecDVEPITCH,vecDVEYAW);
% (old method) eta = LE_vec_local(:,:,2)./2;
vecDVEHVSPN = leVecLocal(:,2)./2;
% 
% % Find Leading Edge Sweep
% % arctan(LE X local component/ LE Y local component)
% (old method) phi_LE = atan(LE_vec_local(:,:,1)./LE_vec_local(:,:,2));
vecDVELESWP = atan(leVecLocal(:,1)./leVecLocal(:,2));
% Find Trailing Edge Sweep
% Project TE Points onto DVE plane
% (TE_Left / TE_Right) (CP)                   (DVE_norm)
% q(x,y,z) TE point | p(a,b,c) Control Point | n(d,e,f) DVE normal
% q_proj = q - dot(q-p,n)*n
% (old method) TE_Left_proj = TE_Left-repmat(dot(TE_Left-CP,DVE_norm,3),1,1,3).*DVE_norm;
teLeftProj = P4 - repmat(dot(P4-matCENTER,matDVENORM,2),1,3).*matDVENORM;
% (old method) TE_Right_proj = TE_Right-repmat(dot(TE_Right-CP,DVE_norm,3),1,1,3).*DVE_norm;
teRightProj = P3 - repmat(dot(P3-matCENTER,matDVENORM,2),1,3).*matDVENORM;
% (old method) TE_vec_proj = TE_Right_proj - TE_Left_proj;
teVecProj = teRightProj-teLeftProj;

% Rotate the Projected TE on DVE to local reference frame
% arctan(Projected TE local X component/Projected TE local Y component)
% (old method) TE_vec_proj_local = fcnGLOBSTAR3D( TE_vec_proj,nu,epsilon,psi );
teVecProjLocal = fcnGLOBSTAR(teVecProj,vecDVEROLL,vecDVEPITCH,vecDVEYAW);
% (old method) phi_TE = atan(TE_vec_proj_local(:,:,1)./TE_vec_proj_local(:,:,2));
vecDVETESWP = atan(teVecProjLocal(:,1)./teVecProjLocal(:,2));


% Compute DVE Mid-chord Sweep
% Average of LE and TE Sweep
% (old method) phi_MID = (phi_LE+phi_TE)./2;
vecDVEMCSWP = (vecDVELESWP+vecDVETESWP)./2;

% Calculating Area
% (old method) Area = eta.*xsi.*4;
vecDVEAREA = vecDVEHVCRD.*vecDVEHVSPN.*4;






% output matDVE, index list which describes the DVE coordinates along with
% vertices location in matVLST
% (old method) verticeList = [LECoordL;LECoordR;TECoordR;TECoordL];
verticeList = [P1;P2;teRightProj;teLeftProj];

%(VAP3)  add vecDVEWING infomation to verticeList before running 'unique'
% this ensures vertices from different vehicle/type(rotor vs wing) 
% will not get WELDED
try
verticeList = [verticeList,double(reshape(repmat(vecDVEWING,1,4),[],1))];
catch
end

[matVLST,idxVLST,matDVE] = unique(verticeList,'rows');

% (VAP3) remove vecDVEWING info from matVLST for output
matVLST = matVLST(:,1:3);
matDVE = reshape(matDVE,valNELE,4);

end

