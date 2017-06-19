function [ matVEHUVW, matVEHROT, vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK )
%FCNINITVEHICLE Summary of this function goes here
%   Inititalize velocites and rotation angles of vehicle

vecVEHYAW = vecVEHTRK + vecVEHBETA;
vecVEHPITCH = vecVEHFPA + vecVEHALPHA;

matVEHROT = deg2rad([vecVEHROLL, vecVEHPITCH, vecVEHYAW]);

[x,y,z] = sph2cart(deg2rad(vecVEHTRK),deg2rad(vecVEHFPA),vecVEHVINF);

matVEHUVW = [-x,y,z];


end

