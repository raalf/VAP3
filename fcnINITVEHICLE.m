function [ matVEHUVW, matVEHROT, matVEHROTRATE, ...
    vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, vecVEHALPHA, ...
    vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK, vecVEHRADIUS )
%FCNINITVEHICLE Summary of this function goes here
%   Inititalize velocites and rotation angles of vehicle

vecVEHYAW = vecVEHTRK + vecVEHBETA;
vecVEHPITCH = vecVEHFPA + vecVEHALPHA;

matVEHROT = deg2rad([vecVEHROLL, vecVEHPITCH, vecVEHYAW]);



% [x,y,z] = sph2cart(deg2rad(vecVEHTRK),deg2rad(vecVEHFPA),vecVEHVINF);


% matVEHUVW = [-x,y,z];

for n = 1:length(vecVEHVINF)
    matVEHUVW(n,:) = -[vecVEHVINF(n) 0 0] * angle2dcm(deg2rad(vecVEHTRK(n)),deg2rad(vecVEHROLL(n)),deg2rad(vecVEHFPA(n)),'ZXY');
end


matVEHROTRATE = zeros(length(vecVEHVINF),3);
matVEHROTRATE(~isnan(vecVEHRADIUS),3) = -matVEHUVW(~isnan(vecVEHRADIUS),1)./vecVEHRADIUS(~isnan(vecVEHRADIUS)); %L/(2*pi*r) = theta/(2*pi)





end

