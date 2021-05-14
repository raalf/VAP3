clear,clc
beam_nodes = [0 0 0; 0 1 0; 0 2 0];
node_norm = [0 0 1; 0 0 1; 0 0 1];
figure(69)
clf
plot3(beam_nodes(:,1),beam_nodes(:,2),beam_nodes(:,3),'-ok','linewidth',1.5,'markersize',8)
box on
hold on

def = [0; 0.2; 0.5];
elevation = [0; asin((def(2:3)-def(1:2))./[1;1])];

[x,y,z] = sph2cart([0; pi/2; pi/2], elevation,[0; 1; 1]);

beam_nodes = beam_nodes + [x,y,z];
beam_nodes(:,3) = cumsum(beam_nodes(:,3));
plot3(beam_nodes(:,1),beam_nodes(:,2),beam_nodes(:,3),'--^r','linewidth',1.5,'markersize',8)