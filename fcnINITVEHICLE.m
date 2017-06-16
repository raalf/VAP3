function [ vecVEHUVW, vecVEHROT ] = fcnINITVEHICLE( vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK )
%FCNINITVEHICLE Summary of this function goes here
%   Inititalize velocites and rotation angles of vehicle

vecVEHYAW = vecVEHTRK + vecVEHBETA;
vecVEHPITCH = vecVEHFPA + vecVEHALPHA;

vecVEHROT = deg2rad([vecVEHROLL, vecVEHPITCH, vecVEHYAW]);

[x,y,z] = sph2cart(deg2rad(vecVEHTRK),rad2deg(vecVEHFPA),vecVEHVINF);

vecVEHUVW = [-x,y,z];


end

