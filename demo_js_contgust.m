V = 40;    % m/s flight speed
dt = 0.0101; % s sampling time
T = 1.25;    % s time length
sigma = 0.25; % m/s turbulence intensity
L = 50;     % m length scale
[cgust,t,Cgust,f] = fcnCONTGUST(V,dt,T,sigma,L);