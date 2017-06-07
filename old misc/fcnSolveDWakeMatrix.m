function [FW] = fcnSolveDWakeMatrix(FW, Aircraft, D_wake, R_wake, Temp,start)

    % Solving of the D-Matrix!!!! San Diego 2016-01-10 by W.J.M.B and T.D.K
    % #FreeWakeAcrossAmerica
    Coeff = D_wake\R_wake;
    
    % Saving the information with all the other DVE info
    % Remember: A,B,C of DVE is equal to 3*DVE Index - 2 to 3*DVE Index
    
for h = start:Temp.timestep+1
    for i = 1:Aircraft.General.Panels
        for j = 1:length(FW.Panels(i).WakeDVE(h).Index)
            indx = FW.Panels(i).WakeDVE(h).Index(j);
%             FW.Panels(i).DVE.A(j,1) = Coeff(2*indx-2);
            
            FW.Panels(i).WakeDVE(h).B(j) = Coeff(2*indx-1);
            FW.Panels(i).WakeDVE(h).C(j) = Coeff(2*indx);
            FW.Panels(i).WakeDVE(h).A(j) = FW.Panels(i).WakeDVE(Temp.timestep+1).K(j) - ((0.5*2/3)* FW.Panels(i).WakeDVE(h).C(j)*FW.Panels(i).WakeDVE(h).eta(j)*FW.Panels(i).WakeDVE(h).eta(j));
            %fprintf('A: %f B: %f C: %f\n',FW.Panels(i).DVE.A(j,1), FW.Panels(i).DVE.B(j,1), FW.Panels(i).DVE.C(j,1));
        end
    end
end
    
end
