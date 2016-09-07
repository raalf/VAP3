function [D] = fcnAssembleWingD(FW, Aircraft)
% 2016-01-08 I think this is all good. Not 100% if it will work for a wing
% defined right-to-left, but I think it should. It currently works for a
% wing built from left-to-right, I think. T.K.

D = zeros(sum([FW.Panels.n])*FW.m*3,sum([FW.Panels.n])*FW.m*3); % Preaccolation for Turbo-Boost in performance

row = 1;

for i = 1:Aircraft.General.Panels
    indx = FW.Panels(i).Edge1; % Index of the chordwise row of elements at the start (Edge 1) of the Panel
    indx_start = FW.Panels(i).Edge1(1)-1; % Subtract this from the element indx so we can call this element from the structure
    
    % Stepping through each column of DVEs from 1 to n for that panel
    % (Again, we are working with the elements in the chordwise direction
    % all at one time for every n
    for ii = 1:FW.Panels(i).n
        % 220 inside the panel elements
        if ii < FW.Panels(i).n; % Only performing this if we aren't at the very outter row of DVEs
            
            for j = 1:size(indx)
                cols(j,:) = 3*indx(j)-2:3*(indx(j)+1); % Matrix containing all of the column numbers in which we will put the various A, B, C values, for DVEs in indx (and indx+1, because this is 220)
            end
            
            for jj = 1:length(cols(:,1)) % Going through the DVEs in indx, which are all the m for this specific n (we are going through the DVEs in the chordwise direction for every n)
                % Applying conditions of constant circulation (A1 + B1eta1 + C1eta1^2 - A2 + B2eta2 - C2eta2^2 = 0)
                % and vorticity (derivative of the previous eqn with respect to eta)
                D(row,cols(jj,:)) = [1 FW.Panels(i).DVE.eta(indx(jj)-indx_start) FW.Panels(i).DVE.eta(indx(jj)-indx_start)^2 -1 FW.Panels(i).DVE.eta(indx(jj)-indx_start+1) -(FW.Panels(i).DVE.eta(indx(jj)-indx_start+1)^2)];
                D(row+1,cols(jj,:)) = [0 1 2*FW.Panels(i).DVE.eta(indx(jj)-indx_start) 0 -1 2*(FW.Panels(i).DVE.eta(indx(jj)-indx_start+1))];
                row = row + 2;
            end
        end
        
        indx = indx+1;
    end
    
    clear cols
    
    % This section enters into the D-matrix the equations for the DVEs on
    % the symmetry plane (vorticity is 0 at the symmetry plane). It does
    % this by first seeing what edge is on the symmetry plane (from the
    % input structure)
    if isempty(Aircraft.Surface(i).Sym) == 0
        if Aircraft.Surface(i).Sym == 2 % If Edge 2 is on the symmetry plane
            indx = FW.Panels(i).Edge2; % Getting index of all panels in this column (1xm), same as above
        else
            indx = FW.Panels(i).Edge1; % If not Edge 2, must be Edge 1
        end
        for j = 1:size(indx)
            cols(j,:) = 3*indx(j)-2:3*indx(j); % Columns corresponding to the DVE elements
            D(row,cols(j,:)) = [0 etaSign(Aircraft.Surface(i).Sym)*1 2*FW.Panels(i).DVE.eta(indx(j)-indx_start)]; % Enter in vorticity equation (same as above)
            row = row + 1;
        end
    end
    clear cols
end

% For every joint in FW.Joints, we will create boundary conditions.
% These joints are either free ends (circulation = 0), or boundaries
% between two panels (circulation, vorticity constant, same as above),
% or boundaries between three panels (circulation, vorticity, curvature
% all constant)
for j = 1:length(FW.Joint)
    clear cols
    clear indx
    
    % Free-end condition
    if length(FW.Joint(j).Panel) == 1
        
        if FW.Joint(j).Edge == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel).Edge1(1)-1;  

        for jj = 1:size(indx)              
            ed = 1; % First panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];    
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
            end
            row = row + 1;
        end
        
        clear cols
        
        % Boundary between two panels
    elseif length(FW.Joint(j).Panel) == 2
        
        % To simplify this, I do it one panel at a time, not updating the
        % row number until both panels have been entered
        
        % Starting with the first panel, finding if its Edge 1 or Edge 2
        if FW.Joint(j).Edge(1) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(1)).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(1)).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel(1)).Edge1(1)-1;
        
        % Entering in the first panel equations, circulation and vorticity
        for jj = 1:size(indx)              
            ed = 1; % First panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 1 2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 1 -2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
            end
            row = row + 2;
        end
        row = row - (jj*2);
                
        clear cols
        
        % Now the second panel
        if FW.Joint(j).Edge(2) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(2)).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(2)).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel(2)).Edge1(1)-1;
        
        % Entering in the second panel equations, circulation and vorticity
        for jj = 1:size(indx)              
            ed = 2; % second panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 -2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
            end
            row = row + 2;
        end
        
        clear cols
        
        % Boundary between three panels
    elseif length(FW.Joint(j).Panel) == 3
        
        % To simplify this, I do it one panel at a time, not updating the
        % row number until both panels have been entered
        
        % Starting with the first panel, finding if its Edge 1 or Edge 2
        if FW.Joint(j).Edge(1) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(1)).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(1)).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel(1)).Edge1(1)-1;
        
        % Entering in the first panel equations, circulation and vorticity and curvature
        for jj = 1:size(indx)              
            ed = 1; % First panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 1 2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 2];
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 1 -2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 2];
            end
            row = row + 3;
        end
        row = row - (jj*3);
        
        clear cols
        
        % Now the second panel
        if FW.Joint(j).Edge(2) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(2)).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(2)).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel(2)).Edge1(1)-1;
        
        % Entering in the second panel equations, circulation and vorticity
        % and curvature
        % NOT updating the row counter
        for jj = 1:size(indx)              
            ed = 2; % First panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 -2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 -2];
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 -2];
            end
            row = row + 3;
        end
        row = row - (jj*3);
        
        clear cols
        
        % Now the third panel
        if FW.Joint(j).Edge(3) == 2 % If the right edge is the one with the condition
            indx = FW.Panels(FW.Joint(j).Panel(3)).Edge2;
        else % Otherwise it must be the left
            indx = FW.Panels(FW.Joint(j).Panel(3)).Edge1;
        end
        
        indx_start = FW.Panels(FW.Joint(j).Panel(3)).Edge1(1)-1;
        
        % Entering in the second panel equations, circulation and vorticity
        % and curvature
        % NOT updating the row counter
        for jj = 1:size(indx)              
            ed = 3; % First panel
            if FW.Joint(j).Edge(ed) == 2 % If the right edge is the one with the condition
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 -2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 -2];
            else % Otherwise it must be the left 
                cols(jj,:) = 3*indx(jj)-2:3*indx(jj); % Columns corresponding to the DVE elements
                D(row,cols(jj,:)) = [-1 FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start) -FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)^2];
                D(row+1,cols(jj,:)) = [0 -1 2*FW.Panels(FW.Joint(j).Panel(ed)).DVE.eta(indx(jj)-indx_start)];
                D(row+2,cols(jj,:)) = [0 0 -2];
            end
            row = row + 3;
        end
        
        clear cols
        
        % Too many panels are joined here, not enough equations!
    else
        disp('Too many panels at joint ', i)
    end
    
end
clear row jj j indx_start indx ii i cols


end

