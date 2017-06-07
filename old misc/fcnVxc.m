function [ Vxc, root_bending, highspeed_cd ] = fcnVxc(Aircraft, Performance, Output, Conditions)

%% Analysis

% section cl * y location * density * 0.5 * section area * V_inf^2
cl_dist = Output.Alpha_5(end).cl;
y_dist = Output.Alpha_5(end).y;
S_dist = Output.Alpha_5(end).S;
vindx = find([Performance.Alpha] == 5);
V_inf = Performance(vindx).Vinf;

root_bending = sum(cl_dist.*y_dist.*Conditions.Density.*0.5.*S_dist.*(V_inf^2));

% Drag coefficient at 100 kts

highspeed_cd = interp1([Performance.Vinf],[Performance.CD],51,'linear','extrap');

% Smoothing data
[LDfit, ~] = createFit([Performance.Alpha], [Performance.LD]);
[CLfit, ~] = createFit([Performance.Alpha], [Performance.CL]);
[CDfit, ~] = createFit([Performance.Alpha], [Performance.CD]);
[Vinffit, ~] = createFit([Performance.Alpha], [Performance.Vinf]);
[Cdifit, ~] = createFit([Performance.Alpha], [Performance.CDi]);

% Vxc

range_vxc = 1.5:0.25:13.5;
CL = CLfit(range_vxc);
CD = CDfit(range_vxc);
LD = LDfit(range_vxc);
Vcruise = Vinffit(range_vxc);
wglide = Vcruise.*(CD./CL);
[~, LDindex] = max(LD);



Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*Aircraft.Reference.Weight/(Aircraft.Reference.S*Conditions.Density);

k = 1;

for wmaxth = 2:0.5:8
    
    j = 1;
    
    for i = LDindex:size(CL)
        wclimb(j,1) = MaxClimb(CL(i), CD(i), Rrecip, wmaxth, WSroh);
        j = j + 1;
    end
    
    [wclimbMAX, indexWC] = max(wclimb);
    
    for i = 1:size(CL)
        V(i,1) = (Vcruise(i)*wclimbMAX)/(wglide(i)+wclimbMAX);
    end
    
    [VxcMAX, cruiseIndex] = max(V);
    invVxcMAX(k,1) = 1/VxcMAX;
    Vxc(k,:) = [wmaxth VxcMAX];
    k = k + 1;
    
end

invVxcMAX_low = invVxcMAX(1,1);
invVxcMAX_med = invVxcMAX(ceil(end/2),1);
invVxcMAX_high = invVxcMAX(end,1);

out = [invVxcMAX_low invVxcMAX_med invVxcMAX_high root_bending highspeed_cd];


end

