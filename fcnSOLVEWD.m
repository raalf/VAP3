function [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN)
% This function solves the D R system of equations and
% returns the coefficients A1 A2 B1 B2 C3 in a matrix, with
% rows corresponding to HDVE numbers

% INPUT:
%   matD - D-matrix, ? x NELE*5 matrix of boundary conditions
%   vecR - Resultant, ? x 1 vector of solutions to the system of equations
% OUTPUT:
%   matCOEFF - NELE x 5 x 1 matrix of (A1, A2, B1, B2, C3) coefficients for each HDVE  

matWCOEFF(:,2:3) = reshape(matWD\vecWR,2,valWNELE,1)';

matWCOEFF(:,1) = vecWKGAM - (1/3).*matWCOEFF(:,3).*(vecWDVEHVSPN.^2);

end

