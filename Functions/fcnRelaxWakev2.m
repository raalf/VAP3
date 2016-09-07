function [FW, Temp] = fcnRelaxWake(FW,Temp,Aircraft)
% 	//relaxing the wake:
% 	//	1. computes local induced velocity at side edges of DVEs
% 	//	2. displaces of side edges of DVEs
% 	//	3. computes new ref. pt of wake DVE
% 	//	4. computes new eta, nu, epsilon, and psi, as well as new xsi

jj = 1;

% I HAVE A CONFESSION TO MAKE. THIS MAY NOT WORK WITH SYMMETRY. I AM
% FOLLOWING ALONG WITH WHAT HE IS DOING IN WAKE_GEOMETRY.CPP IN RELAX_WAKE.
% HERE I AM USING THE SAME NOTATION AS HIM, LEFT AND RIGHT. BUT LETS BE
% HONEST WITH EACH OTHER, THIS WAS NOT WHAT WE HAVE BEEN DOING. WE
% HAVE BEEN USING EDGE 1 AND 2.
% BUT IT IS SUNNY OUTSIDE, AND CHARLIE IS SLEEPING QUIETLY IN HIS BED. YET
% STILL I TOIL. 
% LETS FIX THIS IN POST, WHAT DO YOU SAY, ALTON???
% T.D.K, CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA, 92110. 2016-01-21



for h = Temp.timestep+1:-1:2 % loop through timesteps (for indexing, we start at 1 instead of 0. Starting at 2 to ignore timestep 0)
    jj = 1;
    for i = 1:Aircraft.General.Panels
        % Locally induced velocities at left side of wake DVE
        for j = 1:length(FW.Panels(i).WakeDVE(h).Index)
            P = FW.Panels(i).WakeDVE(h).xleft(j,:);
            Temp.wDVEleft(h).w_ind(jj,:) = fcnDVEInducedVelocity(FW, Aircraft, Temp, P);
            Temp.wDVEleft(h).Panel(jj) = i;
            Temp.wDVEleft(h).sIndex(jj) = j;
            
            %             % w_ind(spanwise location, xyz, timestep)
            %             w_ind(jj,:,h-1) = fcnDVEInducedVelocity(FW, Aircraft, Temp, P);
            jj = jj + 1;
        end
    end
    
    % Locally induced velocities at the free-tips of the wings
    jj = 1;
    for j = 1:length(FW.Joint)
        if length(FW.Joint(j).Panel) == 1
            P = FW.Panels(FW.Joint(j).Panel).WakeDVE(h).xright(end,:);
            % uright(Timestep, Wingtip#, Wing #) Wingtip# may be > 1 for
            % split-tips
            Temp.uright(h).w_ind(jj,:) = fcnDVEInducedVelocity(FW, Aircraft, Temp, P);
            Temp.uright(h).Panel(jj) = FW.Joint(j).Panel;
            Temp.uright(h).sIndex(jj) = jj;
            Temp.uright(h).Wing(jj) = FW.Panels(FW.Joint(j).Panel).Wing;
            jj = jj + 1;
            
        end
    end
end

% -------------------------------------------------------------------------
% Displacing points

% OUR XLEFT, XRIGHT (?) POINTS FROM ABOVE AND BELOW ARE A LITTLE OFF FROM
% BRAMESFELD'S, I SUSPECT BECAUSE I AM USING THE TRUE XLEFT AND XRIGHT, NOT
% THE PLANAR PROJECTIONS FROM ALTONS (REMEMBER WE FOUND THE PROJECTED
% TRAILING EDGE?? YEA YEA??)
% SO I'M JUST GONNA KEEP ON KEEPIN ON AND WE'LL SEE WHAT SHAKES LOOSE.
% #FIXITINPOST
% T.D.K, COURTYARD STUDY LOUNGE, STUDENT CENTER, UNIVERSITY OF CALIFORNIA
% - IRVINE, IRVINE, CALIFORNIA, USA, 92697. 2016-01-23

% T.D.K, RYERSON UNIVERSITY, TORONTO, ONTARIO, CANADA, M5B 2K3. 2016-02-05

% Uncommenting the one will relax row immediately behind TE of wing, do we 
% want this or no?
% T.D.K, RYERSON UNIVERSITY, TORONTO, ONTARIO, CANADA, M5B 2K3. 2016-05-06


for h = 2:Temp.timestep+1
    
    % WHYYYYY  ARE WE NOW INPUTING ALL THE SAME
    % INTO DISPLACE FUNCTION??? WHY, OH WHY GREAT BRAMESFELD???
    if h == 2 % Second oldest row of wakeDVEs (Only current and upstream ??)
        for j = 1:length(Temp.wDVEleft(h).w_ind(:,1))
            P = FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:);
            FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:) = fcnDisplacePoint(FW, Temp.wDVEleft(h).w_ind(j,:), Temp.wDVEleft(h).w_ind(j,:), Temp.wDVEleft(h).w_ind(j,:), P);
%             fprintf('\nxLEFT after 1: %f %f %f', FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:));
        end
        
        % displacing wakeDVE wingtip edges
        for j = 1:length(Temp.uright(h).Panel)
            P = FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:);
            FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:) = fcnDisplacePoint(FW, Temp.uright(h).w_ind(j,:), Temp.uright(h).w_ind(j,:), Temp.uright(h).w_ind(j,:), P);
        end
        
    elseif h == Temp.timestep+1 % Youngest row of wakeDVEs (no upstream elements, only current and downstream)
        
        % WHYYYYY  ARE WE NOW INPUTING [CURRENT, [0 0 0], LAST_TIMESTEP]
        % INTO DISPLACE FUNCTION??? WHY, OH WHY GREAT BRAMESFELD???
        
        % Displacing xleft points for all wakeDVEs
        for j = 1:length(Temp.wDVEleft(h).w_ind(:,1))
            P = FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:);
            %fprintf('xleft before: %f %f %f\n', P(1), P(2), P(3));
            FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:) = fcnDisplacePoint(FW, Temp.wDVEleft(h).w_ind(j,:), [0 0 0], Temp.wDVEleft(h-1).w_ind(j,:), P);
            %fprintf('\nxLEFT after 2: %f %f %f', FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:));
        end
        
        % displacing wakeDVE wingtip edges
        for j = 1:length(Temp.uright(h).Panel)
            P = FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:);
            FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:) = fcnDisplacePoint(FW, Temp.uright(h).w_ind(j,:), [0 0 0], Temp.uright(h-1).w_ind(j,:), P);
        end
    else
        % Displacing xleft points for all wakeDVEs
        for j = 1:length(Temp.wDVEleft(h).w_ind(:,1))
            P = FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:);
            %fprintf('xleft before: %f %f %f\n', P(1), P(2), P(3));
            FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:) = fcnDisplacePoint(FW, Temp.wDVEleft(h-1).w_ind(j,:), Temp.wDVEleft(h).w_ind(j,:), Temp.wDVEleft(h+1).w_ind(j,:), P);
            %fprintf('\nxLEFT after 3: %f %f %f', FW.Panels(Temp.wDVEleft(h).Panel(j)).WakeDVE(h).xleft(Temp.wDVEleft(h).sIndex(j),:));
        end
        
        % displacing wakeDVE wingtip edges
        for j = 1:length(Temp.uright(h).Panel)
            P = FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:);
            FW.Panels(Temp.uright(h).Panel(j)).WakeDVE(h).xright(end,:) = fcnDisplacePoint(FW, Temp.uright(h-1).w_ind(j,:), Temp.uright(h).w_ind(j,:), Temp.uright(h+1).w_ind(j,:), P);
        end
    end
end

% -------------------------------------------------------------------------

% Calculating new reference point xo

for h = Temp.timestep+1:-1:1
    for i = 1:Aircraft.General.Panels
        for j = 1:length(FW.Panels(i).WakeDVE(h).Index)
            
            % Updating the new right components inside of this panel
            if j < length(FW.Panels(i).WakeDVE(h).Index)
                FW.Panels(i).WakeDVE(h).xright(j,:) = FW.Panels(i).WakeDVE(h).xleft(j+1,:);
            elseif j == length(FW.Panels(i).WakeDVE(h).Index)
                % Updating the new right components if between panels
                for mm = 1:length(FW.Joint)
                    if length(FW.Joint(mm).Panel) == 2 % If its a 220 joint
                        [indx1,indx2] = find(FW.Joint(mm).Panel == i);
                        
                        if ~isempty(indx1) && FW.Joint(mm).Edge(indx1) == 2 
                            [indx,~] = find(FW.Joint(mm).Panel ~= i);
                            FW.Panels(i).WakeDVE(h).xright(j,:) = FW.Panels(FW.Joint(mm).Panel(indx)).WakeDVE(h).xleft(1,:);
                        end   
                        
                    elseif length(FW.Joint(mm).Panel) == 3 % If its a split in the wing (333 case)
                        
                        [indx1,indx2] = find(FW.Joint(mm).Panel == i);
                        
                        if ~isempty(indx1) && FW.Joint(mm).Edge(indx1(1)) == 2
                            [indx,~] = find(FW.Joint(mm).Panel ~= i); 
                            FW.Panels(i).WakeDVE(h).xright(j,:) = FW.Panels(FW.Joint(mm).Panel(indx(1))).WakeDVE(h).xleft(1,:);
                        end                          
                        
                    end
                end
            end
            

            
            % This isn't the real norm, it is temporarily being stored here
            [FW.Panels(i).WakeDVE(h).xo(j,:), FW.Panels(i).WakeDVE(h).u(j,:), FW.Panels(i).WakeDVE(h).norm(j,:)] = fcnNewxo(FW, FW.Panels(i).WakeDVE(h).xleft(j,:), FW.Panels(i).WakeDVE(h).xright(j,:), FW.Panels(i).WakeDVE(h).xo(j,:));
            
            if h == Temp.timestep+1
                % nu, eps, psi, xsi from trailing edge of wing
                US_DVE.nu = FW.Panels(i).DVE.roll(j);
                US_DVE.epsilon = FW.Panels(i).DVE.pitch(j);
                US_DVE.psi = FW.Panels(i).DVE.yaw(j);
                US_DVE.xsi = FW.Panels(i).DVE.xsi(j);
                US_DVE.xo = FW.Panels(i).DVE.xo(j,:);
                
                % New roll, pitch, yaw, norm, sweep (LE, MID, TE), xsi
                [FW.Panels(i).WakeDVE(h).norm(j,:), FW.Panels(i).WakeDVE(h).roll(j), FW.Panels(i).WakeDVE(h).pitch(j), FW.Panels(i).WakeDVE(h).yaw(j), ...
                    FW.Panels(i).WakeDVE(h).xsi(j), FW.Panels(i).WakeDVE(h).phiLE(j), FW.Panels(i).WakeDVE(h).phiMID(j), FW.Panels(i).WakeDVE(h).phiTE(j), ...
                    FW.Panels(i).WakeDVE(h).eta(j)] = fcnNewEtaNuEpsPsiXsi(Temp, US_DVE, FW.Panels(i).WakeDVE(h).xo(j,:), FW.Panels(i).WakeDVE(h).norm(j,:), ...
                    FW.Panels(i).WakeDVE(h).phiLE(j), FW.Panels(i).WakeDVE(h).phiMID(j), FW.Panels(i).WakeDVE(h).phiTE(j));
                
            elseif h == 1
                US_DVE.nu = FW.Panels(i).WakeDVE(2).roll(j);
                US_DVE.epsilon = FW.Panels(i).WakeDVE(2).pitch(j);
                US_DVE.psi = FW.Panels(i).WakeDVE(2).yaw(j);
                US_DVE.xsi = FW.Panels(i).WakeDVE(2).xsi(j);
                US_DVE.xo = FW.Panels(i).WakeDVE(2).xo(j,:);
                US_DVE.phiTE = FW.Panels(i).WakeDVE(2).phiTE(j);
                US_DVE.eta = FW.Panels(i).WakeDVE(2).eta(j);
                US_DVE.u = FW.Panels(i).WakeDVE(2).u(j,:);
                
                [FW.Panels(i).WakeDVE(h).norm(j,:), FW.Panels(i).WakeDVE(h).roll(j), FW.Panels(i).WakeDVE(h).pitch(j), FW.Panels(i).WakeDVE(h).yaw(j), ...
                    FW.Panels(i).WakeDVE(h).phiLE(j), FW.Panels(i).WakeDVE(h).phiMID(j), FW.Panels(i).WakeDVE(h).phiTE(j), ...
                    FW.Panels(i).WakeDVE(h).eta(j), FW.Panels(i).WakeDVE(h).xo(j,:), FW.Panels(i).WakeDVE(h).u(j,:)] = fcnNewWakeDVE0(Temp, US_DVE, FW.Panels(i).WakeDVE(1).xsi(j));
                
            else
                % nu, eps, psi, xsi from upstream element
                US_DVE.nu = FW.Panels(i).WakeDVE(h+1).roll(j);
                US_DVE.epsilon = FW.Panels(i).WakeDVE(h+1).pitch(j);
                US_DVE.psi = FW.Panels(i).WakeDVE(h+1).yaw(j);
                US_DVE.xsi = FW.Panels(i).WakeDVE(h+1).xsi(j);
                US_DVE.xo = FW.Panels(i).WakeDVE(h+1).xo(j,:);
                
                % New roll, pitch, yaw, norm, sweep (LE, MID, TE), xsi
                [FW.Panels(i).WakeDVE(h).norm(j,:), FW.Panels(i).WakeDVE(h).roll(j), FW.Panels(i).WakeDVE(h).pitch(j), FW.Panels(i).WakeDVE(h).yaw(j), ...
                    FW.Panels(i).WakeDVE(h).xsi(j), FW.Panels(i).WakeDVE(h).phiLE(j), FW.Panels(i).WakeDVE(h).phiMID(j), FW.Panels(i).WakeDVE(h).phiTE(j), ...
                    FW.Panels(i).WakeDVE(h).eta(j)] = fcnNewEtaNuEpsPsiXsi(Temp, US_DVE, FW.Panels(i).WakeDVE(h).xo(j,:), FW.Panels(i).WakeDVE(h).norm(j,:), ...
                    FW.Panels(i).WakeDVE(h).phiLE(j), FW.Panels(i).WakeDVE(h).phiMID(j), FW.Panels(i).WakeDVE(h).phiTE(j));
            end
            
            clear US_DVE
            
            % Assigning new singfct, 1% of the tip element half-span
            % Using the shorter tip if there are more than one
            % This will probably only work with one split in the wing
            % Or maybe none? I have no clue
            % T.D.K KING ST. E, TORONTO, ONTARIO, CANADA, M5A 1K5
            
            % Bill wrote this loop, I'm just stealing it from generateWakeDVEs
            for k = 1:size(FW.Joint,2)
                if size(FW.Joint(k).Panel) == 1
                    tippanel = FW.Joint(k).Panel;
                    %what wing is this the tip panel of?
                    wing = FW.Panels(tippanel).Wing;
                    %does this match the current panel's wing?
                    if FW.Panels(i).Wing == FW.Panels(tippanel).Wing
                        break                  
                    end
                end
            end

            FW.Panels(i).WakeDVE(h).singfct(j) = 0.01*FW.Panels(tippanel).WakeDVE(h).eta(end);       
            
        end
    end
end




end


