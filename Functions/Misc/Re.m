%% Section Drag
clc
PanelCount = length(FW.Panels);
MaxSpanwise = max([FW.Panels.n]);

% Calculate V_inf (V = sqrt( 2W / rho S CL))
V_inf = sqrt(2*Aircraft.Reference.Weight/(Conditions.Density*Aircraft.Reference.S*CL));

Re = nan(PanelCount,MaxSpanwise);

for i = 1:PanelCount
    n = FW.Panels(i).n;
    m = FW.Panels(i).m;
    
    xsi = FW.Panels(i).DVE.xsi(1:n);
    nu = Conditions.KinematicViscosity;  %nu = mu/rho
   
    % Re = V_inf*2*surfacePtr[m].xsi*tempS*info.m; (From C++) 
    FW.Panels(i).Re = transpose(V_inf*2*xsi*m/nu);
    Re(i,1:n) = FW.Panels(i).Re;
end

clear xsi nu i n m PanelCount MaxSpanwise


maxRe = max(Re(:));
minRe = min(Re(:));








