% Draw DVEs
clear
clc
load('Misc/Cirrus_DVE.mat');

figure(1)
for n = 1:length(FW.Panels)
    DVE = FW.Panels(n).DVE;
    scatter3(DVE.xo(:,1),DVE.xo(:,2),DVE.xo(:,3),'ro');
    hold on
    quiver3(DVE.xo(:,1),DVE.xo(:,2),DVE.xo(:,3),DVE.norm(:,1),DVE.norm(:,2),DVE.norm(:,3));
    
    
end

hold off
axis equal