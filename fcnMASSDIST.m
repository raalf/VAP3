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

SURF.vecWINGIYY = interp1(SURF.vecSTRUCTSPNDIST,INPU.vecJT,SURF.matCENTER(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1,2)).*(2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)); % Wing mass values at each DVE control point y-location

% INPU.vecVEHCG = sum(VEHI.vecWINGMASS(1)*(VEHI.vecWINGCG(1,:)-INPU.matVEHORIG) + VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
%     + VEHI.vecFUSEMASS.*(VEHI.vecFUSECG-INPU.matVEHORIG) + VEHI.vecWINGMASS(2:end).*(VEHI.vecWINGCG(2:end,:)-INPU.matVEHORIG),1)./...
%     (sum(VEHI.vecWINGMASS,1) + VEHI.vecPAYLMASS + VEHI.vecFUSEMASS);

INPU.vecVEHCG = sum(sum(VEHI.vecWINGMASS.*(VEHI.vecWINGCG-INPU.matVEHORIG),1) + VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
    + sum(VEHI.vecFUSEMASS,1).*(VEHI.vecFUSECG-INPU.matVEHORIG),1)./...
    (sum(VEHI.vecWINGMASS,1) + VEHI.vecPAYLMASS + sum(VEHI.vecFUSEMASS,1));

INPU.vecVEHCG = INPU.vecVEHCG + INPU.matVEHORIG;
COND.vecVEHWEIGHT = (VEHI.vecPAYLMASS + sum(VEHI.vecFUSEMASS,1) + sum(VEHI.vecWINGMASS,1))*9.81;

SURF.vecWINGIYY = SURF.vecWINGIYY.*(2*SURF.vecDVEHVSPN(SURF.vecDVELE(SURF.vecWINGTYPE==1)==1)) + sum(SURF.vecVEHMASS.*(SURF.vecWINGCG(:,1)-INPU.vecVEHCG(:,1)).^2 + SURF.vecVEHMASS.*(SURF.vecWINGCG(:,3)-INPU.vecVEHCG(:,3)).^2,1);

VEHI.vecFUSEIYY = sum(VEHI.vecFUSEMASS.*(VEHI.vecFUSEMASSLOC(:,1)-INPU.vecVEHCG(:,1)).^2 + VEHI.vecFUSEMASS.*(VEHI.vecFUSEMASSLOC(:,3)-INPU.vecVEHCG(:,3)).^2,1);

%% Compute vehicle pitch moment of inertia (Iyy) about CG

% Compute vectors between component CG locations and vehicle CG
r_wingcg = fcnGLOBSTAR(VEHI.vecWINGCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecWINGCG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecWINGCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecWINGCG,1),1)));
r_paylcg = fcnGLOBSTAR(VEHI.vecPAYLCG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecPAYLCG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecPAYLCG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecPAYLCG,1),1)));
% r_fusecg = fcnGLOBSTAR(VEHI.vecFUSECG - INPU.vecVEHCG, deg2rad(COND.vecVEHROLL*ones(size(VEHI.vecFUSECG,1),1)), deg2rad(COND.vecVEHALPHA*ones(size(VEHI.vecFUSECG,1),1)), deg2rad(COND.vecVEHBETA*ones(size(VEHI.vecFUSECG,1),1)));

% r_wingcg = SURF.vecWINGCG - INPU.vecVEHCG;
% r_tailcg = VEHI.vecWINGCG(2,:) - INPU.vecVEHCG;
% r_paylcg = VEHI.vecPAYLCG - INPU.vecVEHCG;
% r_fusecg = VEHI.vecFUSECG - INPU.vecVEHCG;
% 
% VEHI.vecIYY = sum(VEHI.vecWINGMASS(1)*((r_wingcg(1,1)).^2+(r_wingcg(1,3)).^2)...
%     + VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
%     + VEHI.vecFUSEMASS.*((r_fusecg(:,1)).^2+(r_fusecg(:,3).^2)) +...
%     VEHI.vecWINGMASS(2:end).*(((r_wingcg(2:end,1)).^2+(r_wingcg(2:end,3)).^2)),1);
% Iyy = zeros(3);
% 
% Right wing
% for i = 1:length(r_wingcg)
%     rr_wing = [r_wingcg(i,2).^2 + r_wingcg(i,3).^2, -r_wingcg(i,1).*r_wingcg(i,2), -r_wingcg(i,1).*r_wingcg(i,2);...
%         -r_wingcg(i,2).*r_wingcg(i,1), r_wingcg(i,1).^2 + r_wingcg(i,3).^2, -r_wingcg(i,2).*r_wingcg(i,3);...
%         -r_wingcg(i,3).*r_wingcg(i,1), -r_wingcg(i,3).*r_wingcg(i,2), r_wingcg(i,1).^2 + r_wingcg(i,2).^2];
%     Iyy = Iyy + SURF.vecVEHMASS(i).*rr_wing;
% end
% 
% Left wing
% r_wingcg(:,2) = -r_wingcg(:,2);
% for i = 1:length(r_wingcg)
%     rr_wing = [r_wingcg(i,2).^2 + r_wingcg(i,3).^2, -r_wingcg(i,1).*r_wingcg(i,2), -r_wingcg(i,1).*r_wingcg(i,2);...
%         -r_wingcg(i,2).*r_wingcg(i,1), r_wingcg(i,1).^2 + r_wingcg(i,3).^2, -r_wingcg(i,2).*r_wingcg(i,3);...
%         -r_wingcg(i,3).*r_wingcg(i,1), -r_wingcg(i,3).*r_wingcg(i,2), r_wingcg(i,1).^2 + r_wingcg(i,2).^2];
%     Iyy = Iyy + SURF.vecVEHMASS(i).*rr_wing;
% end
% 
% rr_tail = [r_tailcg(1,2).^2 + r_tailcg(1,3).^2, -r_tailcg(1,1).*r_tailcg(1,2), -r_tailcg(1,1).*r_tailcg(1,2);...
%     -r_tailcg(1,2).*r_tailcg(1,1), r_tailcg(1,1).^2 + r_tailcg(1,3).^2, -r_tailcg(1,2).*r_tailcg(1,3);...
%     -r_tailcg(1,3).*r_tailcg(1,1), -r_tailcg(1,3).*r_tailcg(1,2), r_tailcg(1,1).^2 + r_tailcg(1,2).^2];
% 
% rr_payl = [r_paylcg(1,2).^2 + r_paylcg(1,3).^2, -r_paylcg(1,1).*r_paylcg(1,2), -r_paylcg(1,1).*r_paylcg(1,2);...
%     -r_paylcg(1,2).*r_paylcg(1,1), r_paylcg(1,1).^2 + r_paylcg(1,3).^2, -r_paylcg(1,2).*r_paylcg(1,3);...
%     -r_paylcg(1,3).*r_paylcg(1,1), -r_paylcg(1,3).*r_paylcg(1,2), r_paylcg(1,1).^2 + r_paylcg(1,2).^2];
% 
% rr_fuse = [r_fusecg(1,2).^2 + r_fusecg(1,3).^2, -r_fusecg(1,1).*r_fusecg(1,2), -r_fusecg(1,1).*r_fusecg(1,2);...
%     -r_fusecg(1,2).*r_fusecg(1,1), r_fusecg(1,1).^2 + r_fusecg(1,3).^2, -r_fusecg(1,2).*r_fusecg(1,3);...
%     -r_fusecg(1,3).*r_fusecg(1,1), -r_fusecg(1,3).*r_fusecg(1,2), r_fusecg(1,1).^2 + r_fusecg(1,2).^2];
% 
% VEHI.matI = Iyy + VEHI.vecWINGMASS(2).*rr_tail + VEHI.vecPAYLMASS.*rr_payl + VEHI.vecFUSEMASS.*rr_fuse;

VEHI.vecIYY = sum(sum(VEHI.vecWINGMASS(2).*((r_wingcg(2,1)).^2+(r_wingcg(2,3)).^2),1)...
    + 2*sum(SURF.vecWINGIYY,1) + VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
    + VEHI.vecFUSEIYY,1);

% VEHI.vecFUSETAILMASS = VEHI.vecPAYLMASS + VEHI.vecFUSEMASS + VEHI.vecWINGMASS(2);
% 
% VEHI.vecFUSETAILCG = sum(VEHI.vecPAYLMASS.*(VEHI.vecPAYLCG-INPU.matVEHORIG)...
%     + VEHI.vecFUSEMASS.*(VEHI.vecFUSECG-INPU.matVEHORIG),1)./...
%     (VEHI.vecWINGMASS(2) + VEHI.vecPAYLMASS + VEHI.vecFUSEMASS);
% 
% VEHI.vecFUSEIYY = sum(VEHI.vecPAYLMASS.*((r_paylcg(:,1)).^2 + (r_paylcg(:,3)).^2) ...
%     + VEHI.vecFUSEMASS.*((r_fusecg(:,1)).^2+(r_fusecg(:,3).^2)),1);

end