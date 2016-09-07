function [Aircraft, Conditions, FW, Output, Performance] = Main_PerfCode2015beta(Aircraft, Conditions, FW)

% clc
% clear

Temp.isverbose = 1;
if isdeployed == 0
    addpath('Functions');
end
warning off

%% Header
if Temp.isverbose == 1
    disp('===========================================================================');
    disp('Versace Aero Performance (Based on FreeWake 2015)');
    disp('Running Version 201x.xx  .                             .');
    disp('Includes stall model    //                             \\');
    disp('No trim solution       //                               \\');
    disp('                      //                                 \\');
    disp('                     //                _._                \\');
    disp('                  .---.              .//|\\.              .---.');
    disp('         ________/ .-. \_________..-~ _.-._ ~-..________ / .-. \_________');
    disp('                 \ ~-~ /   /H-     `-=.___.=-''     -H\   \ ~-~ /');
    disp('                   ~~~    / H          [H]          H \    ~~~');
    disp('                         / _H_         _H_         _H_ \');
    disp('                           UUU         UUU         UUU');
    disp('===========================================================================');
end

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Getting Aircraft info from FW .txt file
% load('Aircraft/Horst-dOne.mat')
% load('Aircraft/StandardCirrus.mat');
% load('Aircraft/designThree.mat');
% load('Aircraft/designFour-Trouble.mat');
% Conditions.Alpha = 8;

% [FW, Aircraft, Conditions] = FWreader('input-LowTherm.txt');
% Aircraft.Surface(1).Sym = 1;
% Conditions.Alpha = 5;
% Aircraft.Surface(6).Sym = 1;

% FW.Relax = 1;
% FW.Maxtime = 10;
% FW.Deltae = 0;
% Aircraft.Surface(3).Sym = 1;


%% Creating temporary values used throughout the program calculations
Temp.DBL_EPS = 1e-14;
Temp.ZERO = 1e-5;
Temp.timestep = -1;
Aircraft.Reference.AR = Aircraft.Reference.b^2/Aircraft.Reference.S;
deltae = 10000;
Performance = [];

%% Creating Output Structure
Output = struct;

%% Figuring out which panels meet with each other
FW = panelBoundary(Aircraft,FW);

%% Splitting Panels into DVEs
% Discretizing each panel of the wing, one panel at a time
for i = 1:Aircraft.General.Panels
    %     [Aircraft, FW] = generateDVEs(Aircraft, FW, i);
    [Aircraft, FW] = generateDVEs_v2(Aircraft, FW, i);
end

% Applying singfct to all surface panels
for i = 1:Aircraft.General.Panels
    [Aircraft, FW] = generateSingfct(Aircraft, FW, i);
end

clear i

%% D Matrix for the Wing BCs
D = fcnAssembleWingD(FW, Aircraft);

%% D Matrix - Adding Kinematic Conditions
D = fcnDVEKinCond(FW, Aircraft, Temp, D);

%% BEGINNING ALPHA LOOP --------------------------------------------------

% % All AoAs to compute
% if Conditions.Alpha(2)-Conditions.Alpha(1) > 0
%     Alpha = Conditions.Alpha(1):Conditions.Alpha(3):Conditions.Alpha(2);
% else
%     Alpha = Conditions.Alpha(1);
% end

Alpha = Conditions.Alpha;

% clear Alpha
% Alpha = [5 7 1];
counter = 1;



for h = 1:length(Alpha)
    
    if Temp.isverbose == 1
        fprintf('\nAlpha = %0.2f \n',Alpha(h));
    end
    
    % Preparing the Output.Alpha-specific Field Name (where we will save
    % all relevant timestep and output information
    s = strcat('Alpha_',strrep(strrep(num2str(Alpha(h)), '.', '_'),'-','N'));
    Output.(s) = [];
    Temp.AI = s; % Alpha Index, used to reference the correct Output.Alhpa file throughout this procedure
    Temp.Alpha = Alpha(h);
    
    % Computing freestream velocity vector
    % This only works for fixed wing aircraft!! Not rotors!!
    Temp.u = [FW.Uinf*cosd(Alpha(h))*cosd(Conditions.Beta) FW.Uinf*sind(Conditions.Beta) FW.Uinf*sind(Alpha(h))*cosd(Conditions.Beta)];
    
    R = fcnDVEResultant(FW, Temp, Aircraft, D);
    
    FW = fcnSolveDMatrix(FW, Aircraft, D, R);
    
    Output = savetimestep(FW,Temp,Output); % Saving the Timestep -1 case
    
    %% BEGINNING TIMESTEP LOOP -------------------------------------------
    
    while (deltae > FW.Deltae || Temp.timestep < FW.Mintime) && Temp.timestep < FW.Maxtime
        
        Temp.timestep = Temp.timestep + 1;
        
        if Temp.isverbose == 1
            fprintf('%d ', Temp.timestep);
        end
        
        FW = fcnMove_Wing(FW, Temp, Aircraft);
        
        %squirt out wake
        for i = 1:Aircraft.General.Panels
            for j = 1:FW.Panels(i).n
                FW = generateWakeDVEs(FW,Output,Temp,i,j);
            end
        end
        
        % Getting new wake vorticity coefficients
        [D_wake, R_wake] = fcnNewWakeVorticity(FW, Aircraft, Temp,Temp.timestep+1);
        [FW] = fcnSolveDWakeMatrix(FW, Aircraft, D_wake, R_wake, Temp,Temp.timestep+1);
        
        % Update surface vorticity, get new resultant
        R = fcnDVEResultant(FW, Temp, Aircraft, D);
        FW = fcnSolveDMatrix(FW, Aircraft, D, R);
        
        
        % Update wake vorticity
        [D_wake, R_wake] = fcnNewWakeVorticity(FW, Aircraft, Temp,1);
        [FW] = fcnSolveDWakeMatrix(FW, Aircraft, D_wake, R_wake, Temp,1);
        
        %Relax wake
        if FW.Relax == 1 && Temp.timestep > 1
            % Update wake vorticity
            [FW, Temp] = fcnRelaxWake(FW,Temp,Aircraft);
            [D_wake, R_wake] = fcnNewWakeVorticity(FW, Aircraft, Temp,1);
            [FW] = fcnSolveDWakeMatrix(FW, Aircraft, D_wake, R_wake, Temp,1);
        end
        
        %             CL, CDi
        Output = savetimestep(FW,Temp,Output);
        N_force = fcnSurfaceDVENormalForces(FW,Temp,Aircraft);
        Output = fcnDVEWingNormalForces(FW, Aircraft, Temp, Conditions, N_force, Output);
        [Output, D_force] = fcnInducedDVEDrag(FW, Aircraft, Temp,Output);
        
        %% cl, cd, cn distribution
        % This hasn't been validated against FW yet (only vaguely, seems good to me)
        count = 1;
        
        for kk = 1:Aircraft.General.Panels
            
            for k = 1:FW.Panels(kk).n
                
                cw_row = FW.Panels(kk).Edge1 + (k-1); % chordwise row of elements
                
                local_row = cw_row - FW.Panels(kk).Edge1(1)+1;
                S(count) = sum(FW.Panels(kk).DVE.area(local_row));
                
                tempS = 2/(FW.Uinf*FW.Uinf*S(count));
                
                cn(count) = sum(N_force(cw_row,5) + N_force(cw_row,6))*tempS;
                cl(count) = sum(N_force(cw_row,1) + N_force(cw_row,2))*tempS;
                cy(count) = sum(N_force(cw_row,3) + N_force(cw_row,4))*tempS;
                cdi(count) = D_force(count)*tempS;
                
                % Works only for m = 1!
                y(count) = FW.Panels(kk).DVE.xo(k,2);
                
                count = count + 1;
            end
            
        end
        
        Output.(Temp.AI)(Temp.timestep+2).cn = cn;
        Output.(Temp.AI)(Temp.timestep+2).cl = cl;
        Output.(Temp.AI)(Temp.timestep+2).cy = cy;
        Output.(Temp.AI)(Temp.timestep+2).cdi = cdi;
        Output.(Temp.AI)(Temp.timestep+2).y = y;
        Output.(Temp.AI)(Temp.timestep+2).S = S;
        
        clear cn cl cy cd tempS S cw_row local_row count
        
        if Temp.timestep > 1
            deltae = Output.(Temp.AI)(Temp.timestep+2).deltae;
        else
            deltae = 10000;
        end
        
        clear N_force
        
    end
    
    Output = savetimestep(FW,Temp,Output);
    deltae = 10000;
    %% Viscous wrapper
    [FW, Output, Performance, Aircraft, Temp] = fcnViscousWrapper(FW, Output, Performance, Aircraft, Temp, Conditions, counter);
    
    % Clearing and resetting for the next run
    Temp.lasttimestep = Temp.timestep;
    Temp.timestep = -1;
    FW.Panels = rmfield(FW.Panels,{'WakeDVE'});
    
    clear cdi D_force k kk
    
    counter = counter + 1;
end % End of alpha loop

% clear h i ii j s Alpha D D_wake R R_wake deltae

%% Plotting
% Coeff = dlmread('AOA5.00.txt','',4,0);

if Temp.isverbose == 1
    disp('Done.');
    plotDVE(Output, Temp, FW);
    
    hFig2 = figure(2);
    clf(2);
    
    hold on
    if FW.m == 1
        for i = 1:Aircraft.General.Panels
            for j = 1:length(FW.Panels(i).DVE.Index)
                
                eta_range = -Output.(Temp.AI)(end).TimestepData(i).DVE.eta(j):Output.(Temp.AI)(end).TimestepData(i).DVE.eta(j)/500:Output.(Temp.AI)(end).TimestepData(i).DVE.eta(j);
                %                 eta_range_plot = FW.Panels(i).DVE.LECoordL(j,2):norm(FW.Panels(i).DVE.LECoordL(j,2)-FW.Panels(i).DVE.LECoordR(j,2))/1000:FW.Panels(i).DVE.LECoordR(j,2);
                eta_range_plot = Output.(Temp.AI)(end).TimestepData(i).DVE.LECoordL(j,2):norm(Output.(Temp.AI)(end).TimestepData(i).DVE.LECoordL(j,2)-Output.(Temp.AI)(end).TimestepData(i).DVE.LECoordR(j,2))/1000:Output.(Temp.AI)(end).TimestepData(i).DVE.LECoordR(j,2);
                
                
                circ = Output.(Temp.AI)(end).TimestepData(i).DVE.A(j) + Output.(Temp.AI)(end).TimestepData(i).DVE.B(j).*eta_range + Output.(Temp.AI)(end).TimestepData(i).DVE.C(j).*(eta_range.^2);
                %                 vort = FW.Panels(i).DVE.B(j) + 2.*FW.Panels(i).DVE.C(j).*(eta_range)
                plot(eta_range_plot, circ,'-b');
                grid on
                
            end
        end
    end
    drawnow
    hold off
end

% toc
% end

