% clc
clear
close all

home_dir = pwd;

% % % % Fixed
% % % % 0.355430335140642,0.0428627646648056,0.0331512545466698,6113.89912090596,0.0117617496618554
% % % % Relaxed
% % % [out, baseline] = fcnBASELINE(home_dir);
% % % 
% % Fixed
% % 0.242090285595947,0.0422745475460615,0.0329607814166055,6153.14012822062,0.0117523495947097
% % Relaxed
% % 0.242208165781066,0.0422731424373225,0.0329602873009932,6153.81964171613,0.0117517829702389
% z = [0.208347 0.148797 0.255929 0.130086 1.00258 0.367333 7.46909 0.723608 0.10309 1.79309 0.43574 7.49991 0.885122 0.0509858 -2.51158 0.575008 7.45888 0.314806 0.0812669 -2.04289 0.585632 7.49963 0.306767 0.0666243 -3.92359];
% [out, split] = fcnOBJECTIVE2(z, home_dir, false);
% 
% % Fixed
% % 0.236978073559861,0.0422285097331260,0.0329382731473526,6151.38082625983,0.0117613675076144
% z = [0.198843 0.359339 0.790319 0.332498 7.46828 0.605641 0.2297 -0.977123 0.38898 7.4995 0.767585 0.0752712 -2.81756];
% [out, single] = fcnOBJECTIVE1(z, home_dir, false);
% save('matlab3.mat');


load('matlab3.mat');
wh = [6 3.5];
%% Vxc
% hFig1 = figure(1);
% clf(1);
% 
% sidx = 5;
% plot(baseline.Vxc(sidx:end,1), 100*(single.Vxc(sidx:end,2) - baseline.Vxc(sidx:end,2))./baseline.Vxc(sidx:end,2), '-ok')
% hold on
% plot(baseline.Vxc(sidx:end,1), 100*(split.Vxc(sidx:end,2) - baseline.Vxc(sidx:end,2))./baseline.Vxc(sidx:end,2), '--b^')
% hold off
% 
% grid minor
% box on
% axis tight
% xlabel('Thermal Core Strength (m/s)','FontSize',10);
% ylabel('Percent Change in V_x_c','FontSize',10)
% legend('Conventional Winglet','Split Winglet','Location','NorthEast')
% saveFig2Latex(hFig1, 'res-vxc.pdf', wh);

%% Sink Rate v Speed
% wh = [3.5 5]
% hFig2 = figure(2);
% clf(2);
% sidx = 0;
% 
% [baseline_sink, ~] = fcnCREATEFIT(baseline.Vcruise.*3.6, baseline.wglide);
% 
% speeds = single.Vinffit(single.range_vxc(1:end - sidx)).*3.6;
% plot(speeds, 100*(single.wglide(1:end - sidx) - baseline_sink(speeds))./baseline_sink(speeds), '-ok')
% hold on
% speeds = split.Vinffit(split.range_vxc(1:end - sidx)).*3.6;
% plot(speeds, 100*(split.wglide(1:end - sidx) - baseline_sink(speeds))./baseline_sink(speeds), '--b^')
% hold off
% 
% % scatter(baseline.Vcruise, baseline.wglide, 'k');
% % hold on
% % scatter(single.Vcruise, single.wglide, 'r');
% % scatter(split.Vcruise, split.wglide, 'b');
% % hold off
% 
% grid minor
% box on
% axis tight
% xlabel('Airspeed (km/h)','FontSize',10);
% ylabel('Percent Change in Sink Rate','FontSize',10)
% legend('Conventional Winglet','Split Winglet','Location','SouthEast')
% saveFig2Latex(hFig2, 'res-wsink.pdf', wh);

%% 2 Sink Rate v Speed
% 
% hFig9 = figure(9);
% clf(9);
% sidx = 0;
% 
% plot(baseline.Vcruise(1:end - sidx).*3.6, -baseline.wglide(1:end - sidx), '-.r')
% hold on
% plot(single.Vcruise(1:end - sidx).*3.6, -single.wglide(1:end - sidx), '-ok')
% plot(split.Vcruise(1:end - sidx).*3.6, -split.wglide(1:end - sidx), '--b^')
% hold off
% % plot(speeds, 100*(split.wglide(1:end - sidx) - baseline_sink(speeds))./baseline_sink(speeds), '--b^')
% 
% 
% % scatter(baseline.Vcruise, baseline.wglide, 'k');
% % hold on
% % scatter(single.Vcruise, single.wglide, 'r');
% % scatter(split.Vcruise, split.wglide, 'b');
% % hold off
% 
% grid minor
% box on
% axis tight
% xlabel('Airspeed (km/h)','FontSize',10);
% ylabel('Sink Rate (m/s)','FontSize',10)
% legend('Base Aircraft','Blended Winglet','Split Winglet','Location','SouthWest')
% saveFig2Latex(hFig9, 'res-wsink2.pdf', wh);

%% CDi v Speed
% wh = [6 3.5]
% hFig3 = figure(3);
% clf(3);
% sidx = 0;
% 
% [baseline_cdi, ~] = fcnCREATEFIT(baseline.Vcruise.*3.6, baseline.Cdifit(baseline.range_vxc));
% 
% speeds = single.Vinffit(single.range_vxc(1:end - sidx)).*3.6;
% plot(speeds, 100*(single.Cdifit(single.range_vxc(1:end - sidx)) - baseline_cdi(speeds))./baseline_cdi(speeds), '-ok')
% hold on
% speeds = split.Vinffit(split.range_vxc(1:end - sidx)).*3.6;
% plot(speeds, 100*(split.Cdifit(split.range_vxc(1:end - sidx)) - baseline_cdi(speeds))./baseline_cdi(speeds), '--b^')
% hold off
% 
% grid minor
% box on
% axis tight
% xlabel('Airspeed (km/h)','FontSize',10);
% ylabel('Percent Change in C_D_i','FontSize',10)
% legend('Conventional Winglet','Split Winglet','Location','SouthWest')
% saveFig2Latex(hFig3, 'res-cdi.pdf', wh);

%% Re
% wh = [6 3.5]
% hFig21 = figure(21);
% clf(21)
% sidx = 1
% 
% idx = single.WING_SWEEP(sidx).vecDVEPANEL > 2;
% nu = 1.46e-05;
% vinf = single.Vcruise;
% len = 2.*single.WING_SWEEP(sidx).vecDVEHVCRD(idx);
% Re = (vinf(sidx).*len)./nu;
% y = single.WING_SWEEP(sidx).WING.vecSPANLOC(idx);
% s = y(1);
% plot(y - s, Re, '-ok')
% 
% idx = split.WING_SWEEP(sidx).vecDVEPANEL > 2 & split.WING_SWEEP(1).vecDVEPANEL < 6;
% nu = 1.46e-05;
% vinf = split.Vcruise;
% len = 2.*split.WING_SWEEP(sidx).vecDVEHVCRD(idx);
% Re = (vinf(sidx).*len)./nu;
% y = split.WING_SWEEP(sidx).WING.vecSPANLOC(idx);
% hold on
% s = y(1);
% plot(y - s, Re, '--^b')
% hold off
% 
% idx = split.WING_SWEEP(sidx).vecDVEPANEL > 5;
% nu = 1.46e-05;
% vinf = split.Vcruise;
% len = 2.*split.WING_SWEEP(sidx).vecDVEHVCRD(idx);
% Re = (vinf(sidx).*len)./nu;
% y = split.WING_SWEEP(sidx).WING.vecSPANLOC(idx);
% hold on
% plot(y - s - (y(1) - s), Re, '-.sb')
% hold off
% 
% grid minor
% box on
% axis tight
% xlabel('Spar Length Along Winglet','Fontsize',10);
% ylabel('Reynolds Number','FontSize',10);
% legend('Conventional Winglet','Split Winglet Upper Surface','Split Winglet Lower Surface','Location','NorthEast')
% saveFig2Latex(hFig21, 'res-re.pdf', wh);

%% CDp v Speed
wh = [3.5 5]
hFig5 = figure(5);
clf(5);
sidx = 0;

[baseline_cdp, ~] = fcnCREATEFIT(baseline.Vcruise.*3.6, baseline.CDfit(baseline.range_vxc) - baseline.Cdifit(baseline.range_vxc));

speeds = single.Vinffit(single.range_vxc(1:end - sidx)).*3.6;
cdp = single.CDfit(single.range_vxc(1:end - sidx)) - single.Cdifit(single.range_vxc(1:end - sidx));
plot(speeds, 100*(cdp - baseline_cdp(speeds))./baseline_cdp(speeds), '-ok')
hold on
speeds = split.Vinffit(split.range_vxc(1:end - sidx)).*3.6;
cdp = split.CDfit(split.range_vxc(1:end - sidx)) - split.Cdifit(split.range_vxc(1:end - sidx));
plot(speeds, 100*(cdp - baseline_cdp(speeds))./baseline_cdp(speeds), '--b^')
hold off

grid minor
box on
axis tight
xlabel('Airspeed (km/h)','FontSize',10);
ylabel('Percent Change in C_D_p','FontSize',10)
legend('Conventional Winglet','Split Winglet','Location','SouthEast')
saveFig2Latex(hFig5, 'res-cdp.pdf', wh);

%% CD v Speed
hFig4 = figure(4);
clf(4);
sidx = 0;

[baseline_cd, ~] = fcnCREATEFIT(baseline.Vcruise.*3.6, baseline.CDfit(baseline.range_vxc));

speeds = single.Vinffit(single.range_vxc(1:end - sidx)).*3.6;
plot(speeds, 100*(single.CD(1:end - sidx) - baseline_cd(speeds))./baseline_cd(speeds), '-ok')
hold on
speeds = split.Vinffit(split.range_vxc(1:end - sidx)).*3.6;
plot(speeds, 100*(split.CD(1:end - sidx) - baseline_cd(speeds))./baseline_cd(speeds), '--b^')
hold off

grid minor
box on
axis tight
xlabel('Airspeed (km/h)','FontSize',15);
ylabel('Percent Change in C_D','FontSize',15)
legend('Conventional Winglet','Split Winglet','Location','SouthEast')
saveFig2Latex(hFig4, 'res-cd.pdf', wh);







