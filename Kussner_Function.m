U = 1;
wg = 0.055;
c = 1;
b = c/2;

t = 0:0.01:50;

s = U*t./b; % Reduced time (s = U*t/b)

kuss = (1 - 0.5*exp(-0.13*s) - 0.5*exp(-s)); % Kussner
kussneralg = (s.^2+s)./(s.^2+2.82*s+0.8);
wag = (1 - 0.165*exp(-0.0455*s) - 0.335*exp(-0.3*s)); % Wagner

Clk = 2*pi*(wg/U)*kuss; % Kussner
Clk2 = 2*pi*(wg/U)*kussneralg;
% Clw = 2*pi*(w/U)*wag; % Wagner

% Cl = 2*pi*(w/U)*(1-0.328*exp(-0.0745.*tau)-0.582*exp(-0.6.*tau));
% Cl = 2*pi*(w/U)*((tau.^2+tau)./(tau.^2+2.82*tau+0.8));

k = 0.18;
C = besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
Sears = C*(besselj(0,k)-1i*besselj(1,k))+1i*besselj(1,k);
w = 2*k*U/c;
f = w/(2*pi);
L = U/f;
Cls = 2*pi*(wg/U)*real(Sears)*real(exp(1i*w*t));

% dk = 1; %step size for integration
% k = [0.00001:dk:100]; %integration variable
% s = [0.001:0.1:25]; %nondim time
% 
% KussnerInt = zeros(1,length(s));
% 
% for sc = 1:length(s)
%     for kc = 1:length(k)        
%         C(kc) = besselh(1,2,k(kc))/(besselh(1,2,k(kc))+1i*besselh(0,2,k(kc)));
%         Sears(kc) = C(kc)*(besselj(0,k(kc))-1i*besselj(1,k(kc)))+1i*besselj(1,k(kc));
%         KussnerInt(sc,kc) = real(Sears(kc)*exp(-1i*k(kc)))*sin(k(kc)*s(sc))/k(kc);
%      end
% end
% 
% Kussner = 2*pi*(w/U)*trapz(k,KussnerInt,2)*2/pi; %trapezoidal integration of sears function for all k's
