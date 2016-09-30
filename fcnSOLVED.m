function [matCOEFF] = fcnSOLVED(matD, vecR, valNELE)
% This function solves the D R system of equations and
% returns the coefficients A1 A2 B1 B2 C3 in a matrix, with
% rows corresponding to HDVE numbers

% INPUT:
%   matD - D-matrix, ? x NELE*5 matrix of boundary conditions
%   vecR - Resultant, ? x 1 vector of solutions to the system of equations
% OUTPUT:
%   matCOEFF - NELE x 5 x 1 matrix of (A1, A2, B1, B2, C3) coefficients for each HDVE  

matCOEFF = matD\vecR;

matCOEFF = reshape(matCOEFF,3,valNELE,1)';

end

