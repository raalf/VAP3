function [Output] = fcnDVEWingNormalForces(FW, Aircraft, Temp, Conditions, N_force, Output)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

q = 0.5*FW.Uinf*FW.Uinf*Aircraft.Reference.S;

Nt_free(1) = 0;
Nt_free(2) = 0;

Nt_ind(1) = 0;
Nt_ind(2) = 0;

Nt_free(1) = sum(N_force(:,1));
Nt_free(2) = sum(N_force(:,3));

Nt_ind(1) = sum(N_force(:,2));
Nt_ind(2) = sum(N_force(:,4));

%fprintf('Nt_free: %f %f\n', Nt_free(1), Nt_free(2));
%fprintf('Nt_ind: %f %f\n', Nt_ind(1), Nt_ind(2));

if FW.Sym == 1 && Conditions.Beta == 0;
    Nt_free(1) = Nt_free(1)*2;
    Nt_ind(1) = Nt_ind(1)*2;
    Nt_ind(2) = Nt_ind(2)*2;
end

CL = (Nt_free(1) + Nt_ind(1))/q;
CLf = Nt_free(1)/q;

CY = (Nt_free(2) + Nt_ind(2))/q;
CYf = Nt_free(2)/q;

CLi = Nt_ind(1)/q;
CYi = Nt_ind(2)/q;

Output.(Temp.AI)(Temp.timestep+2).CL = CL;
Output.(Temp.AI)(Temp.timestep+2).CY = CY;
Output.(Temp.AI)(Temp.timestep+2).CLi = CLi;
Output.(Temp.AI)(Temp.timestep+2).CYi = CYi;
Output.(Temp.AI)(Temp.timestep+2).CLf = CLf;
Output.(Temp.AI)(Temp.timestep+2).CYf = CYf;

end

