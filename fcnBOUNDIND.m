function [aloc, bloc, cloc] = fcnBOUNDIND(endpoints, phi, yaw, fpl)

dbl_eps = 1e-14;

%% Transformation to bound vortex midpoint reference frame
% and other preliminary stuff

eta_translation = mean(endpoints,3);

hspan1 = abs(endpoints(:,1,2) - endpoints(:,1,1))./2;
hspan2 = abs(endpoints(:,2,2) - endpoints(:,2,1))./2;

hspan = hspan1.*cos(yaw) + hspan2.*sin(yaw);

fp_0 = fpl - eta_translation;

eta_0 = fp_0(:,1).*cos(yaw) + fp_0(:,2).*sin(yaw);
xsi_0 = fp_0(:,1).*sin(yaw) + fp_0(:,2).*cos(yaw);
zeta_0 = fp_0(:,3);

%% General case
% From Horstmann's thesis, Appendix 3
a1 = 1 + (tan(phi).^2);
b1 = -(eta_0 + xsi_0.*tan(phi));
c1 = (xsi_0.^2) + (eta_0.^2) + (zeta_0.^2);

re_1 = sqrt((hspan.^2).*a1 - (2.*hspan.*b1) + c1);
re_2 = sqrt((hspan.^2).*a1 + (2.*hspan.*b1) + c1);

G11 = ((a1.*hspan + b1)./((a1.*c1 - b1.^2).*re_2)) - ((a1.*-hspan + b1)./((a1.*c1 - b1.^2).*re_1));
G12 = ((-b1.*hspan - c1)./((a1.*c1 - b1.^2).*re_2)) - ((-b1.*-hspan - c1)./((a1.*c1 - b1.^2).*re_1));
G13 = (((2.*(b1.^2) - a1.*c1)./((a1.*c1 - b1.^2).*a1.*re_2)) + (1./sqrt(a1.^3)).*log(sqrt(a1).*re_2 + a1.*hspan + b1)) - ...
        (((2.*(b1.^2) - a1.*c1)./((a1.*c1 - b1.^2).*a1.*re_1)) + (1./sqrt(a1.^3)).*log(sqrt(a1).*re_1 + a1.*-hspan + b1));
    
% % FreeWake method:
% tempS = ((a1.*c1 - b1.^2).*re_1.*re_2);
% G11	 = 	(a1.*hspan.*(re_1 + re_2) + b1.*(re_1 - re_2))./tempS;
% G12	 = 	(-b1.*hspan.*(re_1 + re_2) - c1.*(re_1 - re_2))./tempS;
% G13	 = 	((2.*b1.*b1 - a1.*c1).*hspan.*(re_1 + re_2) + b1.*c1.*(re_1 - re_2))./(a1.*tempS);
% G13 = 	G13 + log((sqrt(a1).*re_2 + a1.*hspan + b1)./(sqrt(a1).*re_1 - a1.*hspan + b1))./sqrt(a1.*a1.*a1);

a1_xsi = -G11.*zeta_0;
a1_eta = G11.*zeta_0.*tan(phi);
a1_zeta = G11.*(xsi_0 - eta_0.*tan(phi));

b1_xsi = -G12.*zeta_0;
b1_eta = G12.*zeta_0.*tan(phi);
b1_zeta = G12.*(xsi_0 - eta_0.*tan(phi));

c1_xsi = -G13.*zeta_0;
c1_eta = G13.*zeta_0.*tan(phi);
c1_zeta = G13.*(xsi_0 - eta_0.*tan(phi));

%% Special cases
% If the point lies on the bound vortex
idx1 = (abs(xsi_0 - eta_0.*tan(phi)) <= dbl_eps & abs(zeta_0) <= dbl_eps);

len1 = length(idx1);

a1_xsi(idx1) = zeros(len1,1);
a1_eta(idx1) = zeros(len1,1);
a1_zeta(idx1) = zeros(len1,1);

b1_xsi(idx1) = zeros(len1,1);
b1_eta(idx1) = zeros(len1,1);
b1_zeta(idx1) = zeros(len1,1);

c1_xsi(idx1) = zeros(len1,1);
c1_eta(idx1) = zeros(len1,1);
c1_zeta(idx1) = zeros(len1,1);

% If the point lies on the xsi-eta-zeta plane
idx2 = abs(zeta_0) < dbl_eps;

a1_xsi(idx2) = zeros(len1,1);
a1_eta(idx2) = zeros(len1,1);

b1_xsi(idx2) = zeros(len1,1);
b1_eta(idx2) = zeros(len1,1);

c1_xsi(idx2) = zeros(len1,1);
c1_eta(idx2) = zeros(len1,1);

% If the point lies in the plane defined by zeta-axis and bound vortex
idx3 = abs(xsi_0 - eta_0.*tan(phi)) <= dbl_eps;

a1_zeta(idx3) = zeros(len1,1);

b1_zeta(idx3) = zeros(len1,1);

c1_zeta(idx3) = zeros(len1,1);

%% Output format
aloc = [a1_xsi a1_eta a1_zeta];
bloc = [b1_xsi b1_eta b1_zeta];
cloc = [c1_xsi c1_eta c1_zeta];

end







