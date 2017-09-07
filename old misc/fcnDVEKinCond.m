function [D] = fcnDVEKinCond(FW, Aircraft, Temp, D)

% Facsimile of DVE_KinCond from equ_system.cpp

%% Comments from FW:
% 	//This function computes the third part of the D-matrix for DVE's. It
% 	//computes the kinematic condition (tangency requirement at surface DVE)
% 	//at the reference point of a DVE. Specifiaclly, function computes the
% 	//influence coefficients a,b,c that are due to the DVE systems of the
% 	//lifting surface and the wake.  The coefficients multiplied with the
% 	//local bound vortex coefficients, A, B, and C, will yield the locally
% 	//induced velocity. In order to satisfy the kinematic condition, the
% 	//component that is normal to the element surface, cancels the normal
% 	//component of the local velocity that is composed of the free stream
% 	//and the part induced by the wake. This later part determines the
% 	//resutant vector.
% 	//The upper 2/3 of the D-matrix were determined previously in the
% 	//function BoundaryCond and hasn't changed.
% 	//input
% 	//  surfacePtr			geometric information on surface DVE's
% 	//	info				general info
% 	//
% 	//output
% 	//	D					D-matrix, third part, boundary conditions remain
% 	//						the same as in the previous part.


%% Body of Function

sym_condition = FW.Sym; % This symmetry will need to be fixed. Right now, it is either the
row = length(D(:,1))-(sum([FW.Panels.n])*FW.m)+1; % Row number, where to start putting the King Kong conditions

% whole plane is symmetric or nothing. I imagine we want it to be per lifitng surface??

for i = 1:Aircraft.General.Panels
    for j = 1:length(FW.Panels(i).DVE.Index)
        %             cols(j,:) = 3*indx(j)-2:3*(indx(j)+1); % Matrix containing all of the column numbers in which we will put the various A, B, C values, for DVEs in indx (and indx+1, because this is 220)
        
        P = FW.Panels(i).DVE.xo(j,:); % control point of DVE (the point being influenced)
        for ii = 1:Aircraft.General.Panels % looping through influencing elements (all other DVEs)
            for jj = 1:length(FW.Panels(ii).DVE.Index)
                
                DVE_type = 0; % Filament at leading and trailing edge
                singfct = 0;
                [a, b, c] = fcnDVEInfluenceCoeff(FW, Temp, 0, sym_condition, ii, jj, P, singfct, DVE_type,1); % Passing in 0 as the timestep
                cols = 3*FW.Panels(ii).DVE.Index(jj)-2:3*FW.Panels(ii).DVE.Index(jj); % columns for these A, B, C
                
%                  fprintf('DVE_Influence_Coeff A:%f %f %f \n',a(1),a(2),a(3));
%                  fprintf('DVE_Influence_Coeff B:%f %f %f \n',b(1),b(2),b(3));
%                  fprintf('DVE_Influence_Coeff C:%f %f %f \n',c(1),c(2),c(3));
                
                D(row,cols) = [dot(a,FW.Panels(i).DVE.norm(j,:))  dot(b,FW.Panels(i).DVE.norm(j,:)) dot(c, FW.Panels(i).DVE.norm(j,:))];
% D(row,cols) = [dot(a,[0 1 0])  dot(b,[0 1 0]) dot(c, [0 1 0])];
% D(row+1,cols) = [dot(a,[0 0 1])  dot(b,[0 0 1]) dot(c, [0 0 1])];
% D(row+2,cols) = [dot(a,[1 0 0])  dot(b,[1 0 0]) dot(c, [1 0 0])];
                
                clear a b c cols
            end
        end
%         row = row + 3;
        row = row + 1;
    end
end

end

