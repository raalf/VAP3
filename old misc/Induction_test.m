clc
clear


% P - point that is being induced
% xo - reference point of inducing DVE
nu = 0; 
eps = 0;
phiLE = 0;
phiTE = 0;
psi = 0;

xo = [0.5 0.5 0];
P = [0.5 0.5 1];

eta = 0.5;
xsi = 0.5;

DVE_type = 1; % wake DVE, no filament

Temp.DBL_EPS = 1e-14;
singfct = 1;

[a3, b3, c3] = fcnDVEInduction(Temp, P, xo, nu, eps, phiLE, phiTE, psi, eta, xsi, DVE_type, singfct);