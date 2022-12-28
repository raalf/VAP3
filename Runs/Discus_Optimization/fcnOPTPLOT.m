function [] = fcnOPTPLOT(filename)

data = importdata(filename);

% idx = find(data(:,1) >= -0.62);
% data(idx,:) = [];

color = jet(1000);

% EI = [mean(data(:,2:7),2), mean(data(:,8:13),2), mean(data(:,14:20),2)];
% GJ = [mean(data(:,21:26),2), mean(data(:,27:32),2), mean(data(:,33:39),2)];
% ea = [mean(data(:,40:45),2), mean(data(:,46:51),2), mean(data(:,52:58),2)];
% cg = [mean(data(:,59:64),2), mean(data(:,65:70),2), mean(data(:,71:77),2)];
EI = mean(data(:,3:21),2);
GJ = mean(data(:,22:40),2);
ea = mean(data(:,41:59),2);
cg = mean(data(:,60:78),2);
% bending = [data(:,96), data(:,115)]; % [max,min] tip deflections
% torsion = (180/pi).*[data(:,134), data(:,153)]; % [max,min] tip torsion

% Wing is split into 3 different, equal-sized, sections
% Plot results for average of design variables in root section
figure(11)
% subplot(3,2,1)
plot(log10(ea./cg),data(:,1),'ok')
% scatter3(EI(:,1),GJ(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Root Average')
% xlabel('EI (Nm^2)')
% ylabel('GJ (Nm^2)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0 12e5])

% figure(11)
% subplot(3,2,2)
% scatter3(ea(:,1),cg(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Root Average')
% xlabel('EA Location (c)')
% ylabel('CG Location (c)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0.25 0.75])

% Plot results for average of design variables in mid-span section
% figure(12)
% subplot(3,2,3)
% scatter3(EI(:,2),GJ(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Mid-span Average')
% xlabel('EI (Nm^2)')
% ylabel('GJ (Nm^2)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0 12e5])

% figure(13)
% subplot(3,2,4)
% scatter3(ea(:,2),cg(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Mid-span Average')
% xlabel('EA Location (c)')
% ylabel('CG Location (c)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0.25 0.75])

% Plot results for average of design variables in tip section
% figure(14)
% subplot(3,2,5)
% scatter3(EI(:,3),GJ(:,3),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Tip Average')
% xlabel('EI (Nm^2)')
% ylabel('GJ (Nm^2)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0 12e5])

% figure(15)
% subplot(3,2,6)
% scatter3(ea(:,3),cg(:,3),data(:,3),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% title('Tip Average')
% xlabel('EA Location (c)')
% ylabel('CG Location (c)')
% c.Label.String = '-Energy Altitude (m)';
% xlim([0.25 0.75])

% Plot results for maximum bending and torsion deflections at wing tip
% figure(17)
% scatter3(bending(:,1),torsion(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% xlabel('Maximum Tip Deflection (m)')
% ylabel('Maximum Tip Torsion (deg)')
% c.Label.String = '-Energy Altitude (m)';


% figure(16)
% scatter3(ea(:,2),cg(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% xlabel('EA Location (c)')
% ylabel('CG Location (c)')
% c.Label.String = '-Energy Altitude (m)';

