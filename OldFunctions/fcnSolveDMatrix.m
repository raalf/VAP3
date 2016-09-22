function [FW] = fcnSolveDMatrix(FW, Aircraft, D, R)

    % Solving of the D-Matrix!!!! San Diego 2016-01-10 by W.J.M.B and T.D.K
    % #FreeWakeAcrossAmerica
    Coeff = D\R;
    
    % Saving the information with all the other DVE info
    % Remember: A,B,C of DVE is equal to 3*DVE Index - 2 to 3*DVE Index
    for i = 1:Aircraft.General.Panels
        for j = 1:length(FW.Panels(i).DVE.Index)
            indx = FW.Panels(i).DVE.Index(j);
            FW.Panels(i).DVE.A(j) = Coeff(3*indx-2);
            FW.Panels(i).DVE.B(j) = Coeff(3*indx-1);
            FW.Panels(i).DVE.C(j) = Coeff(3*indx);
            
            %fprintf('A: %f B: %f C: %f\n',FW.Panels(i).DVE.A(j), FW.Panels(i).DVE.B(j), FW.Panels(i).DVE.C(j));
        end
    end
    
end

