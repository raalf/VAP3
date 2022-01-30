U = 1;
w = 0.055;
c = 1;
b = c/2;

t = 0:0.01:10;

tau = U*t./b;

Cl = 2*pi*(w/U)*(1 - 0.5*exp(-0.13*tau) - 0.5*exp(-tau));

% Cl = 2*pi*(w/U)*(1-0.328*exp(-0.0745.*tau)-0.582*exp(-0.6.*tau));
% Cl = 2*pi*(w/U)*((tau.^2+tau)./(tau.^2+2.82*tau+0.8));

dk = 0.1; %step size for integration
k = [0.00001:dk:1000]; %integration variable
s = [0.001:0.1:25]; %nondim time

% KussnerInt = zeros(1,length(tau));
% 
% for sc = 1:length(tau)
%     for kc = 1:length(k)        
%         C(kc) = besselh(1,2,k(kc))/(besselh(1,2,k(kc))+1i*besselh(0,2,k(kc)));
%         Sears(kc) = C(kc)*(besselj(0,k(kc))-1i*besselj(1,k(kc)))+1i*besselj(1,k(kc));
%         KussnerInt(sc,kc) = real(Sears(kc)*exp(-1i*k(kc)))*sin(k(kc)*tau(sc))/k(kc);
%      end
% end
% 
% Kussner = 2*pi*(w/U)*trapz(k,KussnerInt,2)*2/pi; %trapezoidal integration of sears function for all k's
