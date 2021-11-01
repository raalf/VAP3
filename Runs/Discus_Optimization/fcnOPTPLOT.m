function [] = fcnOPTPLOT(filename)

data = importdata(filename);

% idx = find(data(:,1) >= -2.37);
% data(idx,:) = [];

color = jet(1000);

EI = [mean(data(:,2:6),2), mean(data(:,7:11),2), mean(data(:,12:65),2)];
GJ = [mean(data(:,17:21),2), mean(data(:,22:26),2), mean(data(:,27:31),2)];
ea = [mean(data(:,32:36),2), mean(data(:,37:41),2), mean(data(:,42:46),2)];
cg = [mean(data(:,47:51),2), mean(data(:,52:56),2), mean(data(:,57:61),2)];
bending = [data(:,76), data(:,91)]; % [max,min] tip deflections
torsion = (180/pi).*[data(:,106), data(:,121)]; % [max,min] tip torsion

% Wing is split into 3 different, equal-sized, sections
% Plot results for average of design variables in root section
figure(10)
subplot(3,2,1)
scatter3(EI(:,1),GJ(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Root Average')
xlabel('EI (Nm^2)')
ylabel('GJ (Nm^2)')
c.Label.String = '-Energy Altitude (m)';
xlim([0 12e5])

% figure(11)
subplot(3,2,2)
scatter3(ea(:,1),cg(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Root Average')
xlabel('EA Location (c)')
ylabel('CG Location (c)')
c.Label.String = '-Energy Altitude (m)';
xlim([0.25 0.75])

% Plot results for average of design variables in mid-span section
% figure(12)
subplot(3,2,3)
scatter3(EI(:,2),GJ(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Mid-span Average')
xlabel('EI (Nm^2)')
ylabel('GJ (Nm^2)')
c.Label.String = '-Energy Altitude (m)';
xlim([0 12e5])

% figure(13)
subplot(3,2,4)
scatter3(ea(:,2),cg(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Mid-span Average')
xlabel('EA Location (c)')
ylabel('CG Location (c)')
c.Label.String = '-Energy Altitude (m)';
xlim([0.25 0.75])

% Plot results for average of design variables in tip section
% figure(14)
subplot(3,2,5)
scatter3(EI(:,3),GJ(:,3),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Tip Average')
xlabel('EI (Nm^2)')
ylabel('GJ (Nm^2)')
c.Label.String = '-Energy Altitude (m)';
xlim([0 12e5])

% figure(15)
subplot(3,2,6)
scatter3(ea(:,3),cg(:,3),data(:,3),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
title('Tip Average')
xlabel('EA Location (c)')
ylabel('CG Location (c)')
c.Label.String = '-Energy Altitude (m)';
xlim([0.25 0.75])

% Plot results for maximum bending and torsion deflections at wing tip
figure(16)
scatter3(bending(:,1),torsion(:,1),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
colormap(color)
c = colorbar;
box on
view(0,90)
xlabel('Maximum Tip Deflection (m)')
ylabel('Maximum Tip Torsion (deg)')
c.Label.String = '-Energy Altitude (m)';


% figure(16)
% scatter3(ea(:,2),cg(:,2),data(:,1),50*ones(length(data),1),data(:,1),'filled','markeredgecolor','k')
% colormap(color)
% c = colorbar;
% box on
% view(0,90)
% xlabel('EA Location (c)')
% ylabel('CG Location (c)')
% c.Label.String = '-Energy Altitude (m)';

