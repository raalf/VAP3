function [INPU, SURF, VEHI, COND] = fcnMASSDIST(INPU, VEHI, SURF, COND)
% Function to compute new CG location after vehicle movement/deformation

VEHI.vecWINGCG(1,:) = sum(SURF.vecVEHMASS.*(SURF.vecWINGCG-INPU.matVEHORIG),1)./(sum(SURF.vecVEHMASS,1)); % Wing CG location relative to vehicle origin

VEHI.vecWINGCG(1,:) = VEHI.vecWINGCG(1,:) + INPU.matVEHORIG;

VEHI.vecWINGMASS(1) = sum(SURF.vecVEHMASS,1)*2;

if any(INPU.vecSYM == 1)
    VEHI.vecWINGCG(:,2) = 0; % Set y component of wing CG to zero for symmetry
end

% INPU.vecVEHCG = sum(VEHI.vecWINGMASS(1)*(VEHI.vecWINGCG(1,:)-INPU.matVEHORIG) + VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
%     + VEHI.vecFUSEMASS.*(VEHI.vecFUSECG-INPU.matVEHORIG) + VEHI.vecWINGMASS(2:end).*(VEHI.vecWINGCG(2:end,:)-INPU.matVEHORIG),1)./...
%     (sum(VEHI.vecWINGMASS,1) + VEHI.vecPAYLMASS + VEHI.vecFUSEMASS);

INPU.vecVEHCG = sum(sum(VEHI.vecWINGMASS.*(VEHI.vecWINGCG-INPU.matVEHORIG),1) + VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
    + VEHI.vecFUSEMASS.*(VEHI.vecFUSECG-INPU.matVEHORIG),1)./...
    (sum(VEHI.vecWINGMASS,1) + VEHI.vecPAYLMASS + VEHI.vecFUSEMASS);

INPU.vecVEHCG = INPU.vecVEHCG + INPU.matVEHORIG;
COND.vecVEHWEIGHT = (VEHI.vecPAYLMASS + VEHI.vecFUSEMASS + sum(VEHI.vecWINGMASS,1))*9.81;

%% Compute vehicle pitch moment of inertia (Iyy) about CG

% Compute vectors between component CG locations and vehicle CG
r_wingcg = fcnGLOBSTAR(VEHI.vecWINGCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecWINGCG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecWINGCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecWINGCG,1),1)));
r_paylcg = fcnGLOBSTAR(VEHI.vecPAYLCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecPAYLCG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecPAYLCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecPAYLCG,1),1)));
r_fusecg = fcnGLOBSTAR(VEHI.vecFUSECG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecFUSECG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecFUSECG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecFUSECG,1),1)));

% VEHI.vecIYY = sum(VEHI.vecWINGMASS(1)*((r_wingcg(1,1)).^2+(r_wingcg(1,3)).^2)...
%     + VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
%     + VEHI.vecFUSEMASS.*((r_fusecg(:,1)).^2+(r_fusecg(:,3).^2)) +...
%     VEHI.vecWINGMASS(2:end).*(((r_wingcg(2:end,1)).^2+(r_wingcg(2:end,3)).^2)),1);

VEHI.vecIYY = sum(sum(VEHI.vecWINGMASS.*((r_wingcg(:,1)).^2+(r_wingcg(:,3)).^2),1)...
    + VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
    + VEHI.vecFUSEMASS.*((r_fusecg(:,1)).^2+(r_fusecg(:,3).^2)),1);

VEHI.vecFUSETAILCG = sum(VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
    + VEHI.vecFUSEMASS.*(VEHI.vecFUSECG-INPU.matVEHORIG),1)./...
    (VEHI.vecWINGMASS(2) + VEHI.vecPAYLMASS + VEHI.vecFUSEMASS);

VEHI.vecFUSEIYY = sum(VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
    + VEHI.vecFUSEMASS.*((r_fusecg(:,1)).^2+(r_fusecg(:,3).^2)),1);

end