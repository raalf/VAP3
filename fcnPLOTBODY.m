clc
clear

VAP_MAIN;

hFig2 = figure(2);
clf(2);

patch('Faces',matDVE,'Vertices',matVLST,'FaceColor','r')







hold on

patch('Faces',matDVE,'Vertices',[matVLST(:,1),-matVLST(:,2),matVLST(:,3)],'FaceColor','r')
quiver3(vecDVECTLPT(:,1), vecDVECTLPT(:,2), vecDVECTLPT(:,3), vecDVENORM(:,1), vecDVENORM(:,2), vecDVENORM(:,3))

hold off

alpha(0.5)

box on
grid on
axis equal

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);