function [qmil_filename, collective_for_camber] = fcnQMILCREATE(temp_name, airfoil_data, blades, thrust, vinf, rpm, diam)
hub_radius = (diam*0.1875)/2;

vref = (2*pi*((diam/2)*0.75))*(rpm/60);
reref = (vref*0.15)/(1.4e-5);
reexp = -0.5;

%% Finding airfoil data at this reference Re
pol = airfoil_data.pol;
Cl  = reshape(pol(:,2,:),[],1);
Cdp = reshape(pol(:,3,:),[],1);
Re  = reshape(pol(:,8,:),[],1);
Alpha = reshape(repmat(airfoil_data.Alpha_Range',1,1,length(airfoil_data.Re_range)),[],1);

idxNans = isnan(Cl) | isnan(Cdp) | isnan(Re);
Cl = Cl(~idxNans);
Cdp = Cdp(~idxNans);
Re = Re(~idxNans);
Alpha = Alpha(~idxNans);

% Lift-related airfoil data for QMIL
Re_lift = scatteredInterpolant(Re,deg2rad(Alpha),Cl,'linear');
alphas = deg2rad([2:8]');
cla = mean(diff(Re_lift(repmat(reref,length(alphas),1), alphas))./diff(alphas));
cl0 = Re_lift(reref,0);

% y = mx + b
% 0 = cl_a.*alpha + cl0
% cl0./cl_a
collective_for_camber = rad2deg(cl0/cla);  

alphas = deg2rad([-5:15]');
[clmax,idx_max] = max(Re_lift(repmat(reref,length(alphas),1), alphas));
[clmin,idx_min] = min(Re_lift(repmat(reref,length(alphas),1), alphas));

% Lift-related airfoil data for QMIL
Re_drag = scatteredInterpolant(Re,deg2rad(Alpha),Cdp,'linear');
[cd0,idx_0] = min(Re_drag(repmat(reref,length(alphas),1), alphas));
clcd0 = Re_lift(reref, alphas(idx_0));

cd2u = (Re_drag(reref, alphas(idx_max)) - cd0)/((clmax - clcd0).^2);
cd2l = (Re_drag(reref, alphas(idx_min)) - cd0)/((clmin - clcd0).^2);

%%
xides = [0 0.5 1];
cldes = [0.6 0.5 0.4];

str_1 = sprintf('temp prop\n\n%d !Nblades\n\n%.4f %.4f ! CL0 CL_A\n%.4f %.4f ! CLmin CLmax\n\n%.4f %.4f %.4f %.4f ! CD0 CD2u CD2l CLCD0\n%.4f %.4f ! Reref REexp\n\n%.4f %.4f %.4f ! XIdes (r over R locations where design cl is specified)\n%.4f %.4f %.4f ! CLdes (specified cl)\n\n',...
    blades, cl0, cla, clmin, clmax, cd0, cd2u, cd2l, clcd0, reref, reexp, xides(1), xides(2), xides(3), cldes(1), cldes(2), cldes(3));

str_2 = sprintf('%.4f ! hub radius(m)\n%.4f ! tip radius (m)\n%.4f ! speed (m/s)\n%.4f ! rpm\n\n', hub_radius, diam/2, vinf, rpm);
str_3 = sprintf('%.4f ! Thrust(N)\n%.4f ! Power(W)\n\n%.4f %.4f ! Ldes KQdes\n\n%d ! Nout number of output stations (optional)', thrust, 0, 0, 0, 30);

str = [str_1, str_2, str_3];

qmil_filename = ['aux_files/', temp_name, '.prop'];

fp = fopen(qmil_filename, 'w');
fprintf(fp, str);
fclose(fp);

end