function [w_ind] = fcnDVEInducedVelocity(FW, Aircraft, Temp, P)
%FCNDVEINDUCEDVELOCITY Summary of this function goes here
%   Detailed explanation goes here

w_surface = fcnSurface_DVE_Vel_Induction ( FW,Aircraft,Temp,P);
w_wake = fcnWake_DVE_Vel_Induction( FW,Aircraft,Temp,P );

w_ind = w_wake + w_surface;

end

