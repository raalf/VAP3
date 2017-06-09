function [D_wake, R_wake] = fcnNewWakeVorticity(FW, Aircraft, Temp,start)
% Updating the vorticity coefficients for the wake, to account for
% stretching. Follows New_Vorticity_Coefficients from wake_geometry.cpp

% This will loosely follow fcnAssembleWingD, the difference being now
% we are in the wake (so m = 1 essentially)

%start is the timestep to start the wake at. This input will either be the current timestep, or it will be the number 1 (to count through all steps. Ask Bill...)
D_wake = zeros(sum([FW.Panels.n])*2,sum([FW.Panels.n])*2); % Preaccolation for Turbo-Boost in performance
R_wake = zeros(sum([FW.Panels.n])*2,1);

row = 1;

for h = start:Temp.timestep+1 % going through all timesteps
    
    for i = 1:Aircraft.General.Panels % going through all panels
        
        for ii = 1:FW.Panels(i).n % going through all spanwise wakeDVEs per panel
            
            indx = FW.Panels(i).WakeDVE(h).Index(ii);
            % 220 inside the panel elements
            if ii < FW.Panels(i).n; % Only performing this if we aren't at the very outter row of DVEs
                
                cols = 2*indx-1:2*(indx+1); % Matrix containing all of the column numbers in which we will put the various A, B, C values, for wakeDVEs in indx (and indx+1, because this is 220)
                
                D_wake(row,cols) = [FW.Panels(i).WakeDVE(h).eta(ii) (2/3)*FW.Panels(i).WakeDVE(h).eta(ii)^2 FW.Panels(i).WakeDVE(h).eta(ii+1) -(2/3)*FW.Panels(i).WakeDVE(h).eta(ii+1)^2];
                D_wake(row+1,cols) = [1 2*FW.Panels(i).WakeDVE(h).eta(ii) -1 2*(FW.Panels(i).WakeDVE(h).eta(ii+1))];
                
                R_wake(row) = -FW.Panels(i).WakeDVE(Temp.timestep+1).K(ii) + FW.Panels(i).WakeDVE(Temp.timestep+1).K(ii+1); % Change in circulation between panels (signs may be wrong??)
                
                R_wake(row+1) = 0; 
                row = row + 2;
            end
        end
        
        clear cols
        
        % This section enters into the D-matrix the equations for the wakeDVEs on
        % the symmetry plane (vorticity is 0 at the symmetry plane). It does
        % this by first seeing what edge is on the symmetry plane (from the
        % input structure)
        if isempty(Aircraft.Surface(i).Sym) == 0
            if Aircraft.Surface(i).Sym == 2 % If Edge 2 is on the symmetry plane
                indx = FW.Panels(i).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
                cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
                D_wake(row,cols) = [etaSign(Aircraft.Surface(i).Sym)*1 2*FW.Panels(i).WakeDVE(h).eta(end)]; % Enter in vorticity equation (same as above)
            else
                indx = FW.Panels(i).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
                cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
                D_wake(row,cols) = [etaSign(Aircraft.Surface(i).Sym)*1 2*FW.Panels(i).WakeDVE(h).eta(1)]; % Enter in vorticity equation (same as above)
            end
            R_wake(row) = 0;
            row = row + 1;
        end
        
        clear cols
        
    end
    



for j = 1:length(FW.Joint)
    clear cols
    clear indx
    
    % Free-end condition
    if length(FW.Joint(j).Panel) == 1
 
        if FW.Joint(j).Edge == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel).WakeDVE(h).eta(end) (2/3)*FW.Panels(FW.Joint(j).Panel).WakeDVE(h).eta(end)^2];
            R_wake(row) = -FW.Panels(FW.Joint(j).Panel).WakeDVE(Temp.timestep+1 ).K(end);
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel).WakeDVE(h).eta(1) (2/3)*FW.Panels(FW.Joint(j).Panel).WakeDVE(h).eta(1)^2];
            R_wake(row) = FW.Panels(FW.Joint(j).Panel).WakeDVE(Temp.timestep+1 ).K(1);
        end

        row = row + 1;
        
        % Boundary between two panels
    elseif length(FW.Joint(j).Panel) == 2
        
        % To simplify this, I do it one panel at a time, not updating the
        % row number until both panels have been entered
        
        % First Panel
        ed = 1;
        if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end) (2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)^2];
            D_wake(row+1,cols) = [1 2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)];
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(end);% unsure
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1) (2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)^2];
            D_wake(row+1,cols) = [1 -2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)];
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(1); %unsure
        end
        
        % Second Panel
        ed = 2;
        if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
            % THIS BLOCK IS TOTALLY SUSPECT, DONT TRUST IT AT ALL. THE
            % RESULTANT STUFF IS SUPER SKETCHY
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)^2];
            D_wake(row+1,cols) = [-1 -2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)];
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(end) - R_wake(row); %not sure about this
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)^2];
            D_wake(row+1,cols) = [-1 2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)];
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(1) - R_wake(row); % not sure about this
        end
        
        R_wake(row+1) = 0;
        row = row + 2;
        
        clear cols
        
        
        
        % Boundary between three panels
    elseif length(FW.Joint(j).Panel) == 3
        
        % To simplify this, I do it one panel at a time, not updating the
        % row number until both panels have been entered
        
        % First Panel
        ed = 1;
        if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end) (2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)^2];
            D_wake(row+1,cols) = [1 2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)];
            D_wake(row+2,cols) = [0 2]; % IS THIS CORRECT?????????????????????????
            R_wake(row) = -FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(end);% unsure
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1) (2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)^2];
            D_wake(row+1,cols) = [1 -2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)];
            D_wake(row+2,cols) = [0 2]; % IS THIS CORRECT?????????????????????????
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(1);% unsure
        end
        
        % Second Panel
        ed = 2;
        if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)^2];
            D_wake(row+1,cols) = [-1 -2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)];
            D_wake(row+2,cols) = [0 -2]; % IS THIS CORRECT?????????????????????????
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(end) - R_wake(row); % not sure about this
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)^2];
            D_wake(row+1,cols) = [-1 2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)];
            D_wake(row+2,cols) = [0 -2]; % IS THIS CORRECT?????????????????????????
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(1) + R_wake(row); % not sure about this
        end
        
        
        % Third Panel
        ed = 3;
        if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(end); % Getting index of all panels in this column (1xm), same as above
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [-FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)^2];
            D_wake(row+1,cols) = [-1 -2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(end)];
            D_wake(row+2,cols) = [0 -2]; % IS THIS CORRECT?????????????????????????
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(end) - R_wake(row); % not sure about this
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).Index(1); % If not Edge 2, must be Edge 1
            cols = 2*indx-1:2*indx; % Columns corresponding to the DVE elements
            D_wake(row,cols) = [FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1) -(2/3)*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)^2];
            D_wake(row+1,cols) = [-1 2*FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(h).eta(1)];
            D_wake(row+2,cols) = [0 -2]; % IS THIS CORRECT?????????????????????????
            
            % THIS WAS FIXED, THANK GOD BILL BISSONNETTE IS A GENIUS.
            % OUR RESULTS NOW MATCH HORSTMANN'S
            % AERODYNAMICS #1 LAB VIBRATIONS LAB, KERR HALL EAST 33
            % RYERSON UNIVERSITY, TORONTO, ONTARIO, CANADA M5B 2K3
            % 2016-05-02
            R_wake(row) = FW.Panels(FW.Joint(j).Panel(ed)).WakeDVE(Temp.timestep+1 ).K(1) + R_wake(row); % not sure about this
        end
        
        R_wake(row+1) = 0;
        R_wake(row+2) = 0;
        row = row + 3;
        
        clear cols
   
        % Too many panels are joined here, not enough equations!
    else
        disp('Too many panels at joint ', i)
    end
    
end

end

clear row jj j indx_start indx ii i cols
end

