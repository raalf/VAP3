function [FW, Output, Performance, Aircraft, Temp] = fcnViscousWrapper(FW, Output, Performance, Aircraft, Temp, Conditions, counter)

Dprof = 0;

q_inf = Aircraft.Reference.Weight/(Output.(Temp.AI)(end).CL*Aircraft.Reference.S);
V_inf = sqrt(2*q_inf/Conditions.Density);

Di = Output.(Temp.AI)(end).CDi*Aircraft.Reference.S*q_inf;

tempS = 1/Conditions.KinematicViscosity;

%% Wing/Horizontal Stabilizer drag
% Should I use cn or cl???
count = 1;
for i = 1:Aircraft.General.Panels
    D_panelprofile(i) = 0;
    for j = 1:FW.Panels(i).n
        Re = V_inf*2*FW.Panels(i).DVE.xsi(j)*tempS*FW.m;
        airfoil = FW.Panels(i).Airfoil;
        % load airfoil data
        Temp.Airfoil = dlmread(strcat('airfoils/airfoil',num2str(airfoil),'.dat'),'\t', 1, 0);
        
        % limits of the Re in our airfoil file
        HiRe = Temp.Airfoil(end,4);
        LoRe = Temp.Airfoil(1,4);
        
        cl = Output.(Temp.AI)(end).cn(count);
        
        % finding section cl_max based on airfoil data
        if Re > HiRe
            if Temp.isverbose == 1
                disp('Re higher than airfoil Re data')
            end
            Re1 = airfoil(end,4);
            temp_var = airfoil(airfoil(:,4)==Re2,2);
            cl_max = temp_var(end);
        elseif Re < LoRe
            if Temp.isverbose == 1
                disp('Re lower than airfoil Re data');
            end
            Re2 = Temp.Airfoil(1,4);
            temp_var = Temp.Airfoil(Temp.Airfoil(:,4)==Re2,2);
            cl_max = temp_var(end);
        else
            Re1 = Temp.Airfoil(Temp.Airfoil(:,4)<Re,4);
            Re1 = Re1(end);
            cl_max1 = Temp.Airfoil(Temp.Airfoil(:,4)<Re,2);
            cl_max1 = cl_max1(end);
            
            temp_var = Temp.Airfoil(Temp.Airfoil(:,4)>Re,4);
            Re2 = temp_var(1);
            temp_var = Temp.Airfoil(Temp.Airfoil(:,4)==(temp_var(1)),2);
            cl_max2 = temp_var(end);
            
            cl_max = interp1([Re1 Re2],[cl_max1 cl_max2], Re);
        end
        
        % correcting the section cl if we are above cl_max
        if cl > cl_max
            if Temp.isverbose == 1
                fprintf('Stall of Panel %d Section %d, cl = %f Re = %0.0f\n',i,j,cl,Re)
            end
            Output.(Temp.AI)(end).cn(count) = 0.825*cl_max; % setting the stalled 2d cl
        end
        
        % determining the drag coefficient corresponding to the
        % uncorrected section cl
        % MATLAB:
%        F = scatteredInterpolant(Temp.Airfoil(:,4), Temp.Airfoil(:,2), Temp.Airfoil(:,3),'nearest'); 
        %cd = F(Re, cl);
      % Octave:
        cd = griddata(Temp.Airfoil(:,4), Temp.Airfoil(:,2), Temp.Airfoil(:,3), Re, cl, 'nearest');
        D_panelprofile(i) = D_panelprofile(i) + (cd*FW.Panels(i).DVE.area(j)*FW.m);
        %fprintf('\nRe: %0.0f\ncl = %f\ncorrected cl = %f\ncd = %f\n', Re, cl, Output.(Temp.AI)(end).cn(count), cd)
        
        count = count + 1;
    end
    
    
    if Aircraft.Surface(i).Sym == 1
        Dprof = Dprof + 2*(D_panelprofile(i)*q_inf);
    else
        Dprof = Dprof + (D_panelprofile(i)*q_inf);
    end
    
end
%% Vertical tail drag
Dvt = 0;
for ii = 1:FW.VerticalStabilizer.Panels
    Re = V_inf*FW.VerticalStabilizer.Geometry(ii,2)*tempS;
    airfoil = FW.VerticalStabilizer.Geometry(ii,4);
    
    % load airfoil data
    Temp.Airfoil = dlmread(strcat('airfoils/airfoil',num2str(airfoil),'.dat'),'\t', 1, 0);
    
    % determining the drag coefficient corresponding to lift
    % coefficient of 0

            % MATLAB:
%    F = scatteredInterpolant(Temp.Airfoil(:,4), Temp.Airfoil(:,2), Temp.Airfoil(:,3),'nearest');
%        cdvt = F(Re, 0);
      % Octave:
        cdvt = griddata(Temp.Airfoil(:,4), Temp.Airfoil(:,2), Temp.Airfoil(:,3), Re, 0, 'nearest');

    Dvt = Dvt + cdvt*FW.VerticalStabilizer.Geometry(ii,3);
end

Dvt = Dvt*q_inf;

%% Fuselage drag

Dfuselage = 0;

tempSS = V_inf*FW.Fuselage.Width*tempS;

for ii = 1:FW.Fuselage.Sections
    Re_fus = (ii-0.5)*tempSS;
    if ii < FW.Fuselage.Turbulent
        cdf = 0.664/sqrt(Re_fus); % Laminar
    else
        cdf = 0.0576/(Re_fus^0.2); % Turbulent
    end
    
    Dfuselage = Dfuselage + cdf*FW.Fuselage.Geometry(ii,2)*pi*FW.Fuselage.Width;
end

Dfuselage = Dfuselage*q_inf;

%% Total Drag
Dmisc = 0;

D = Di + Dprof + Dvt + Dfuselage + Dmisc;

Dint = D*(FW.Interference/100);

D = D + Dint;

%% Adjusting CL for stall
% Should I use cn or cl???
tempS = 0;
count = 1;
for i = 1:Aircraft.General.Panels
    for j = 1:FW.Panels(i).n
        tempS = tempS + Output.(Temp.AI)(Temp.timestep+2).cn(count)*FW.Panels(i).DVE.area(j)*FW.m*cosd(FW.Panels(i).DVE.roll(j));
        count = count + 1;
    end
end

tempS = tempS/Aircraft.Reference.S*2;

if Temp.isverbose == 1
    fprintf('Stall-corrected CL = %f, old was %f\n', tempS, Output.(Temp.AI)(end).CL)
    fprintf('\nCL: %f CDi: %f e: %f\n', Output.(Temp.AI)(end).CL, Output.(Temp.AI)(end).CDi, Output.(Temp.AI)(end).e);
end

Performance(counter).Alpha = Temp.Alpha;
Performance(counter).CLi = Output.(Temp.AI)(end).CL;
Performance(counter).CL = tempS;
Performance(counter).CDi = Output.(Temp.AI)(end).CDi;
Performance(counter).CD = D/(q_inf*Aircraft.Reference.S);
Performance(counter).Vinf = V_inf;
Performance(counter).LD = Performance(counter).CL/Performance(counter).CD;
Performance(counter).Di = Di;
Performance(counter).D = D;
Performance(counter).Dprof = Dprof;
Performance(counter).Dvt = Dvt;
Performance(counter).Dfuselage = Dfuselage;
Performance(counter).Dmisc = Dmisc;


end

