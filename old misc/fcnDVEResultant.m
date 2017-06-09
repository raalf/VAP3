function [R] = fcnDVEResultant(FW, Temp, Aircraft, D)
% Building resultant Vector
% Top 2/3 will be 0 (corresponding to panel BCs)
% Lower 1/3 will be resultants from the King Kong conditions

R = zeros(length(D(:,1)),1); %Preallocating for Turbo-Boost in performance
row = length(D(:,1))-(sum([FW.Panels.n])*FW.m)+1; % Row number, where to start putting the King Kong conditions in the resultant
% row = 21;
for i = 1:Aircraft.General.Panels
    for j = 1:length(FW.Panels(i).DVE.Index)

        if Temp.timestep < 0
            w_extern(1) = Temp.u(1);
            w_extern(2) = Temp.u(2);
            w_extern(3) = Temp.u(3);
        else
            % this is where we find the influence of the wake
            P = FW.Panels(i).DVE.xo(j,:);
            w_wake = fcnWake_DVE_Vel_Induction( FW,Aircraft,Temp,P );
%             fprintf('w_wake: %f %f %f\n', w_wake(1), w_wake(2), w_wake(3));
            w_extern = Temp.u + w_wake;
        end

        R(row) = 4*pi*dot(w_extern, FW.Panels(i).DVE.norm(j,:));
        
%         R(row) = 4*pi*dot(w_extern, [0 1 0]);
%         R(row+1) = 4*pi*dot(w_extern, [0 0 1]);
%         R(row+2) = 4*pi*dot(w_extern, [1 0 0]);

%         fprintf('w_extern: %f %f %f\n', w_extern(1), w_extern(2), w_extern(3));
%         fprintf('Normal: %f %f %f\n', FW.Panels(i).DVE.norm(j,1), FW.Panels(i).DVE.norm(j,2),FW.Panels(i).DVE.norm(j,3));
%         fprintf('Resultant: %f\n', R(row));

% row = row + 3;
        row = row + 1;
    end
end

%% Trying to vectorize, fcnWake_DVE_Vel_Induction needs to be vectorized,
% right now w_extern works for timestep < 0 ( I think). Just vectorize
% induction and it should be good
% for i = 1:Aircraft.General.Panels
%     
%     if Temp.timestep < 0
%         w_extern = repmat(Temp.u,FW.Panels(i).n*FW.Panels(i).m,1);     
%     else
%         w_wake = fcnWake_DVE_Vel_Induction(FW,Aircraft,Temp,FW.Panels(i).DVE.xo);
%         w_extern = Temp.u + w_wake(1:end,:);
%     end
%         
%     R(row:row+2) = 4*pi*dot(w_extern, FW.Panels(i).DVE.norm,2);
%   
%     row = row + 3;
% end

clear i j s

end

