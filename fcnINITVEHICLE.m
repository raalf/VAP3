function [ matVEHUVW, matVEHROT, matVEHROTRATE, ...
    matCIRORIG] = fcnINITVEHICLE( vecVEHVINF, matVEHORIG, vecVEHALPHA, ...
    vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS )
%FCNINITVEHICLE Summary of this function goes here
%   Inititalize velocites and rotation angles of vehicle

vecVEHYAW = vecVEHTRK + vecVEHBETA;
vecVEHPITCH = vecVEHFPA + vecVEHALPHA(end);

matVEHROT = deg2rad([vecVEHROLL, vecVEHPITCH, vecVEHYAW]);

% [x,y,z] = sph2cart(deg2rad(vecVEHTRK),deg2rad(vecVEHFPA),vecVEHVINF);


% matVEHUVW = [-x,y,z];
% velocities UVW of each vehicle
for n = 1:length(vecVEHVINF)
    matVEHUVW(n,:) = -[vecVEHVINF(n) 0 0] * angle2dcm(deg2rad(vecVEHTRK(n)),deg2rad(vecVEHROLL(n)),deg2rad(vecVEHFPA(n)),'ZXY');

% cirling origin of vehicle flight path
    matCIRORIG(n,:) = matVEHORIG(n,:) + [vecVEHRADIUS(n) 0 0] * angle2dcm(deg2rad(vecVEHTRK(n))+pi/2,0,0,'ZXY');

end

% calculate the angular rate of each vehicle if it is doing circling
% flight based on velocity U
matVEHROTRATE = zeros(length(vecVEHVINF),3);
matVEHROTRATE(~isnan(vecVEHRADIUS),3) = -matVEHUVW(~isnan(vecVEHRADIUS),1)./vecVEHRADIUS(~isnan(vecVEHRADIUS)); %L/(2*pi*r) = theta/(2*pi)





end
