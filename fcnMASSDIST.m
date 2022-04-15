function [INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND)
% Function to compute new CG location after vehicle movement/deformation

VEHI.vecWINGCG(1,:) = sum(SURF.vecVEHMASS.*(SURF.vecWINGCG-INPU.matVEHORIG),1)./(sum(SURF.vecVEHMASS,1)); % Wing CG location relative to vehicle origin

VEHI.vecFUSECG = sum(VEHI.vecFUSEMASS.*(VEHI.vecFUSEMASSLOC-INPU.matVEHORIG),1)./(sum(VEHI.vecFUSEMASS,1)); % Fuse CG location relative to vehicle origin

VEHI.vecWINGCG(1,:) = VEHI.vecWINGCG(1,:) + INPU.matVEHORIG;
VEHI.vecFUSECG(1,:) = VEHI.vecFUSECG(1,:) + INPU.matVEHORIG;

VEHI.vecWINGMASS(1) = sum(SURF.vecVEHMASS,1)*2;

if any(INPU.vecSYM == 1)
    VEHI.vecWINGCG(:,2) = 0; % Set y component of wing CG to zero for symmetry
end

SURF.vecWINGIYY = interp1(SURF.vecSTRUCTSPNDIST,INPU.vecJT,SURF.matCENTER(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1,2)).*...
    (2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)); % Wing inertia values at each DVE control point y-location

INPU.vecVEHCG = sum(sum(VEHI.vecWINGMASS.*(VEHI.vecWINGCG-INPU.matVEHORIG),1) + VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
    + sum(VEHI.vecFUSEMASS,1).*(VEHI.vecFUSECG-INPU.matVEHORIG),1)./...
    (sum(VEHI.vecWINGMASS,1) + VEHI.vecPAYLMASS + sum(VEHI.vecFUSEMASS,1));

INPU.vecVEHCG = INPU.vecVEHCG + INPU.matVEHORIG;
COND.vecVEHWEIGHT = (VEHI.vecPAYLMASS + sum(VEHI.vecFUSEMASS,1) + sum(VEHI.vecWINGMASS,1))*9.81; % Total aircraft weight

SURF.vecWINGIYY = SURF.vecWINGIYY + SURF.vecVEHMASS.*(SURF.vecWINGCG(:,1)-INPU.vecVEHCG(:,1)).^2 + SURF.vecVEHMASS.*(SURF.vecWINGCG(:,3)-INPU.vecVEHCG(:,3)).^2; % Parallel-axis theorem for wing inertia about aircraft CG

VEHI.vecFUSEIYY = sum(sum(VEHI.vecFUSEMASS,1).*(VEHI.vecFUSECG(:,1)-INPU.vecVEHCG(:,1)).^2 + sum(VEHI.vecFUSEMASS,1).*(VEHI.vecFUSECG(:,3)-INPU.vecVEHCG(:,3)).^2,1); % Parallel-axis theorem for fuselage inertia about aircraft CG

%% Compute vehicle pitch moment of inertia (Iyy) about CG

% Compute vectors between component CG locations and vehicle CG
r_wingcg = fcnGLOBSTAR(VEHI.vecWINGCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecWINGCG,1),1)),...
    deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecWINGCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecWINGCG,1),1)));
r_paylcg = fcnGLOBSTAR(VEHI.vecPAYLCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecPAYLCG,1),1)),...
    deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecPAYLCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecPAYLCG,1),1)));

VEHI.vecIYY = sum(sum(VEHI.vecWINGMASS(2).*((r_wingcg(2,1)).^2+(r_wingcg(2,3)).^2),1)...
    + 2*sum(SURF.vecWINGIYY,1) + VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
    + VEHI.vecFUSEIYY,1);

end