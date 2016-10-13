function [aloc, bloc, cloc] = fcnVSIND(hspan, phi, fp_0, k)
% This function finds the influence of a semi-infinite vortex sheet on a
% point

% INPUT:
%   hspan - half-span of the DVE
%   phi - sweep of the vortex sheet leading edge
%   fp_0 - vector to the field point from the mid-point of the leading edge
%           of the vortex sheet (origin), rotated into the reference frame of the DVE. 

% OUTPUT:
%   aloc, bloc, cloc - influence coefficients of the vortex sheet on the point

% T.D.K 2016-09-28 ROTHWELL STREET, AURORA, ONTARIO, CANADA, L4G-0V8

dbl_eps = 1e-14;

eta_0 = fp_0(:,2);
xsi_0 = fp_0(:,1);
zeta_0 = fp_0(:,3);

zeta_0sq = zeta_0.*zeta_0;

tanphi = tan(phi);

len = length(eta_0(:,1));

le_vect = xsi_0 - eta_0.*tanphi;

% Eqn A2-12
a2 = 1 + (tanphi.^2);
b2 = le_vect.*tanphi;
c2 = le_vect.^2 + zeta_0sq;
t1 = eta_0 + hspan;
t2 = eta_0 - hspan;

t1s = t1.*t1;
t2s = t2.*t2;

rt_1 = sqrt((t1s).*a2 + 2.*t1.*b2 + c2);
rt_2 = sqrt((t2s).*a2 + 2.*t2.*b2 + c2);

% Eqn A2-5
eps = (le_vect.^2) - (zeta_0sq).*(tanphi).^2;
rho = sqrt(eps.^2 + 4.*(zeta_0sq).*(b2.^2));
beta1 = -sqrt((rho + eps)./2);
beta2 = -sqrt((rho - eps)./2);

% Corrections to beta for special conditions
beta1(0.5.*(rho + eps) <= dbl_eps) = 0;
beta2(0.5.*(rho - eps) <= dbl_eps) = 0;

zetab2 = zeta_0.*b2;
idx_B2 = (abs(zetab2) > dbl_eps);
beta2(idx_B2) = beta2(idx_B2).*(zetab2(idx_B2))./abs(zetab2(idx_B2));

% Eqn A2-8
mu3_1 = a2.*t1 + b2 + sqrt(a2).*rt_1;
mu3_2 = a2.*t2 + b2 + sqrt(a2).*rt_2;

% gamma1 = (1./rho).*(a2.*beta2.*zeta_0 + b2.*beta1);
% gamma2 = (1./rho).*(a2.*beta1.*zeta_0 - b2.*beta2);
% delta1 = (1./rho).*(b2.*beta2.*zeta_0 + c2.*beta1);
% delta2 = (1./rho).*(b2.*beta1.*zeta_0 - c2.*beta2);
% mu1_1 = ((gamma1.*t1 + delta1 - rt_1).^2 + (gamma2.*t1 + delta2).^2)./(k + t1.^2 + zeta_0.^2);
% mu1_2 = ((gamma1.*t2 + delta1 - rt_2).^2 + (gamma2.*t2 + delta2).^2)./(k + t2.^2 + zeta_0.^2);
% mu2_1 = atan(zeta_0./t1) + atan((gamma2.*t1 + delta2)./(gamma1.*t1 + delta1 - rt_1));
% mu2_2 = atan(zeta_0./t2) + atan((gamma2.*t2 + delta2)./(gamma1.*t2 + delta1 - rt_2));

% He changed the above equations to the below ones
gamma1 = (a2.*beta2.*zeta_0 + b2.*beta1);
gamma2 = (a2.*beta1.*zeta_0 - b2.*beta2);
delta1 = (b2.*beta2.*zeta_0 + c2.*beta1);
delta2 = (b2.*beta1.*zeta_0 - c2.*beta2);
mu1_1 = ((gamma1.*t1 + delta1 - rt_1.*rho).^2 + (gamma2.*t1 + delta2).^2);
mu1_2 = ((gamma1.*t2 + delta1 - rt_2.*rho).^2 + (gamma2.*t2 + delta2).^2);
mu2_1 = atan(zeta_0./t1) + atan((gamma2.*t1 + delta2)./(gamma1.*t1 + delta1 - rt_1.*rho));
mu2_2 = atan(zeta_0./t2) + atan((gamma2.*t2 + delta2)./(gamma1.*t2 + delta1 - rt_2.*rho));

%% Implementing special conditions

idx_m21 = (t1s < dbl_eps);
mu2_1(idx_m21) = (pi/2).*abs(zeta_0(idx_m21))./zeta_0(idx_m21) + ...
    atan((gamma2(idx_m21).*t1(idx_m21) + delta2(idx_m21))./ ...
    (gamma1(idx_m21).*t1(idx_m21) + delta1(idx_m21) - rt_1(idx_m21).*rho(idx_m21)));

mu2_1(zeta_0 > 0 & t1 < 0 & t1s > dbl_eps) = mu2_1(zeta_0 > 0 & t1 < 0 & t1s > dbl_eps) + pi;
mu2_1(zeta_0 < 0 & t1 < 0 & t1s > dbl_eps) = mu2_1(zeta_0 < 0 & t1 < 0 & t1s > dbl_eps) - pi;
mu2_1(gamma1.*t1 + delta1 - rt_1.*rho < 0) = mu2_1(gamma1.*t1 + delta1 - rt_1.*rho < 0) + pi;
mu2_1(gamma2.*t1 + delta2 < 0 & gamma1.*t1 + delta1 - rt_1.*rho > 0) = ...
    mu2_1(gamma2.*t1 + delta2 < 0 & gamma1.*t1 + delta1 - rt_1.*rho > 0) + 2*pi;

idx_m22 = (t2s < dbl_eps);
mu2_2(idx_m22) = (pi/2).*abs(zeta_0(idx_m22))./zeta_0(idx_m22) + ...
    atan((gamma2(idx_m22).*t2(idx_m22) + delta2(idx_m22))./ ...
    (gamma1(idx_m22).*t2(idx_m22) + delta1(idx_m22) - rt_2(idx_m22).*rho(idx_m22)));

mu2_2(zeta_0 > 0 & t2 < 0 & ~idx_m22) = mu2_2(zeta_0 > 0 & t2 < 0 & ~idx_m22) + pi;
mu2_2(zeta_0 < 0 & t2 < 0 & ~idx_m22) = mu2_2(zeta_0 < 0 & t2 < 0 & ~idx_m22) - pi;
mu2_2(gamma1.*t2 + delta1 - rt_2.*rho < 0) = mu2_2(gamma1.*t2 + delta1 - rt_2.*rho < 0) + pi;
mu2_2(gamma2.*t2 + delta2 < 0 & gamma1.*t2 + delta1 - rt_2.*rho > 0) = ...
    mu2_2(gamma2.*t2 + delta2 < 0 & gamma1.*t2 + delta1 - rt_2.*rho > 0) + 2*pi;

% Using half span instead of half chord - unconfirmed if this is OK
idx31 = abs(phi) > dbl_eps;
mu3_1(idx31) = 0.0001.*hspan(idx31) + mu3_1(idx31);
mu3_2(idx31) = 0.0001.*hspan(idx31) + mu3_2(idx31);

%%
% 
% G25b = zeros(len,1);
% G25c = zeros(len,1);
% G26a = zeros(len,1);
G21b = zeros(len,1);
G21c = zeros(len,1);

G25b = -0.5.*log((k + zeta_0sq + t2s)./(k + zeta_0sq + t1s));
G25c = -hspan.*log((k + zeta_0sq + t1s).*(k + zeta_0sq + t2s));


idx70 = abs(t1) > dbl_eps;
% G25c(idx70) = G25c(idx70) + t1(idx70).*log(zeta_0(idx70) + t1s(idx70)); before speed boost
G25c70 = G25c + t1.*log(zeta_0 + t1s); %speed boost without index
G25c(idx70) = G25c70(idx70); %speed boost index

idx71 = abs(t2) > dbl_eps;
% G25c(idx71) = G25c(idx71) - t2(idx71).*log(zeta_0(idx71) + t2s(idx71)); before speed boost
G25c71 = G25c - t2.*log(zeta_0 + t2s); %speed boost without index
G25c(idx71) = G25c71(idx71); %speed boost index




% Eqn A2-9
% G25 = (0.5.*log(k + t2.^2 + zeta_0sq)) - (0.5.*log(k + t1.^2 + zeta_0sq));
% G25 = (0.5.*log(t2s + zeta_0sq)) - (0.5.*log(t1s + zeta_0sq)); before speed boost
G25 = (log(t2s + zeta_0sq) - log(t1s + zeta_0sq))./2;



% Eqn A2-3
% G21 = ((beta1./(2.*rho)).*log(mu1_2) + (beta2./rho).*mu2_2) - ((beta1./(2.*rho)).*log(mu1_1) + (beta2./rho).*mu2_1);
G21 = (beta1.*(0.5.*log(mu1_2./mu1_1) - G25) + beta2.*(mu2_2 - mu2_1))./rho;
 
% Eqn A2-4
% G22 = ((1./xsi_0).*(-(beta2./(2.*rho)).*log(mu1_2) + (beta1./rho).*mu2_2)) ...
%     - ((1./xsi_0).*(-(beta2./(2.*rho)).*log(mu1_1) + (beta1./rho).*mu2_1));
G22 = (-beta2.*(0.5.*log(mu1_2./mu1_1) - G25) + beta1.*(mu2_2 - mu2_1))./(rho.*zeta_0);

% Eqn A2-6
lmu3_2 = log(mu3_2);
lmu3_1 = log(mu3_1);
a2cus = sqrt(a2.*a2.*a2);

G23 = ((1./a2).*rt_2 - (b2./a2cus).*lmu3_2) - ((1./a2).*rt_1 - (b2./a2cus).*lmu3_1);

% Eqn A2-7
G24 = ((1./sqrt(a2)).*lmu3_2) - ((1./sqrt(a2)).*lmu3_1);

% Eqn A2-10
G26 = ((1./zeta_0).*atan(t2./zeta_0)) - (1./zeta_0).*atan(t1./zeta_0);
% G26 = atan((t1.*t2 - t1.*zeta_0)./(zeta_0.^2 + t1.*t2)); % Divide by zeta????????????????????????????????????

X2Y2 = (t1.*t2)./(zeta_0sq);
X2 = t2./zeta_0;

idx50 = X2Y2 < -1 & X2 > 0;
G26(idx50) = G26(idx50) + pi;

idx51 = X2Y2 < -1 & X2 < 0;
G26(idx51) = G26(idx51) - pi;

% Eqn A2-11
G27 = t2 - t1;

% Eqn A2-13
b21 = -le_vect;
b22 = (zeta_0sq).*tanphi;
b23 = zeros(len,1);
b24 = -tanphi;
b25 = -ones(len,1);
b26 = zeros(len,1);
b27 = zeros(len,1);

c21 = -2.*((zeta_0sq).*tanphi + eta_0.*le_vect);
c22 = -2.*(zeta_0sq).*(xsi_0 - 2.*eta_0.*tanphi);
c23 = 2.*tanphi;
c24 = 2.*(xsi_0 - 2.*eta_0.*tanphi);
c25 = -2.*eta_0;
c26 = -2.*(zeta_0sq);
c27 = repmat(2,len,1);

% Point is in plane of vortex sheet, but not on bound vortex
idx30 = abs(zeta_0) <= dbl_eps & abs(le_vect) > dbl_eps;
G21(idx30) = 0;
G21b(idx30) = b21(idx30).*beta1(idx30).*(0.5.*log(mu1_2(idx30)./mu1_1(idx30)) + G25b(idx30))./rho(idx30);
G21c(idx30) = b21(idx30).*beta1(idx30).*(eta_0(idx30).*log((mu1_2(idx30)+k(idx30))./(mu1_1(idx30)+k(idx30))) + G25c(idx30))./rho(idx30);
G22(idx30) = 0;
G26(idx30) = 0;

G24(logical(abs(mu3_2 <= dbl_eps | abs(mu3_1) <= dbl_eps))) = 0;

% Point is in plane spanned by bound vortex and zeta axis (only without sweep)
idx20 = abs(xsi_0) <= dbl_eps & abs(zeta_0) > dbl_eps & abs(phi) <= dbl_eps;
% G26a(idx20) = atan(((t2(idx20) - t1(idx20)).*zeta_0(idx20)./(zeta_0(idx20).^2 + t1(idx20).*t2(idx20))));
G21(idx20) = 0;
G22(idx20) = 0;

% idx21 = (idx20 & (zeta_0.^2 + t1.*t2) < 0 & t2./zeta_0 > 0 );
% G26a(idx21) = G26(idx21) + pi;
% 
% idx22 = (idx20 & (zeta_0.^2 + t1.*t2) < 0 & t2./zeta_0 < 0 );
% G26a(idx22) = G26(idx22) - pi;
% 
% G26(idx20) = G26a(idx20)./zeta_0(idx20);

idx23 = (abs(le_vect) <= dbl_eps & abs(zeta_0) <= dbl_eps);
G26(idx23) = 0;

% Eqn A2-2
b2_eta = -zeta_0.*(G21.*b24 + G22.*b21 + G26.*b25);
b2_zeta = G21.*b21 + G22.*b22 + G23.*b23 + G24.*b24 + G25.*b25 + G26.*b26 + G27.*b27;

c2_eta = -zeta_0.*(G21.*c24 + G22.*c21 + G24.*c23 + G25.*c27 + G26.*c25);
c2_zeta = G21.*c21 + G22.*c22 + G23.*c23 + G24.*c24 + G25.*c25 + G26.*c26 + G27.*c27;

b2_xsi = zeros(size(b2_eta));
c2_xsi = zeros(size(c2_eta));

% If point lies on eta-xsi plane
idx40 = abs(zeta_0) <= dbl_eps;
b2_eta(idx40) = 0;
c2_eta(idx40) = 0;
b2_zeta(idx40) = G21b(idx40) + G24(idx40).*b24(idx40) + G25b(idx40);
c2_zeta(idx40) = G21c(idx40) + G23(idx40).*c23(idx40) + G24(idx40).*c24(idx40) + G25c(idx40) + G27(idx40).*c27(idx40);

% If he point falls on a swept leading edge inside the bounds of a sheet
idx60 = abs(zeta_0) <= dbl_eps & abs(le_vect) <= dbl_eps & abs(tanphi) > dbl_eps & abs(hspan) - abs(eta_0) >= -dbl_eps;
b2_zeta(idx60) = 0;
% c2_zeta(idx60) = 0;

% a, b, c in local ref frame
bloc = [b2_xsi b2_eta b2_zeta];
cloc = [c2_xsi c2_eta c2_zeta];
aloc = zeros(size(bloc));

% % If the point lies on a swept leading edge
% idx_LE = abs(zeta_0) <= dbl_eps & abs(xsi_0 - eta_0.*tan(phi)) <= dbl_eps; %& abs(phi) <= dbl_eps;
% bl(idx_LE,:) = zeros(size(bl(idx_LE,:)));
% cl(idx_LE,:) = zeros(size(cl(idx_LE,:)));
% 
% clear idx_LE 

% If the point lies on an unswept leading edge
% a23ind.f - Line 604
idx_LE = abs(zeta_0) <= dbl_eps & abs(xsi_0.^2) <= dbl_eps & abs(phi) <= dbl_eps;
bloc(idx_LE,1:2) = zeros(size(bloc(idx_LE,1:2)));
cloc(idx_LE,1:2) = zeros(size(cloc(idx_LE,1:2)));
% Horstmanns:
% bl(idx_LE,3) = -log((t2(idx_LE) + k(idx_LE))./(t1(idx_LE) + k(idx_LE)));
% cl(idx_LE,3) = -(2.*eta_0(idx_LE).*bl(idx_LE,3) - 2.*(t2(idx_LE) - t1(idx_LE)));

% Reverted to GB's method, seems to work better (less singularities) T.D.K 2016-10-05
% Bramesfelds:
bloc(idx_LE,3) = 0.5.*log((t1s(idx_LE) + k(idx_LE))./(t2s(idx_LE) + k(idx_LE)));
cloc(idx_LE,3) = -4.*hspan(idx_LE) + eta_0(idx_LE).*2.*bloc(idx_LE,3);

end

