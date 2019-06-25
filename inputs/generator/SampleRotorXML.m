% Vogi 5x4.5 rotor
% Using the APC Thin Electric 9x4.5 rotor as a sample
clear,clc

rotor_radius = 0.2286/2; % Radius (m)
airfoil_name = 'MH-117_new';
output_filename = 'APC_E9x4.5';
%% Read illinois geom data
filename = 'APC_E9x4.5_geom.txt';
rotor_data = dlmread(filename,'',1,0);
r_R = rotor_data(:,1);
c_R = rotor_data(:,2);
beta = rotor_data(:,3);

%% Leading edge distribution (x values for each r/R)
% Create from the grabit.m file on a top down photo of the rotor
% le = (-c_R(1)*rotor_radius/2)*(ones(size(r_R,1),1));
tempLE = [0.0903762043168347,0.192658193700100;0.0974874953284367,0.231801282583870;0.101651051582755,0.272655924945316;0.106597062782724,0.325644507408664;0.109433461960709,0.369630040278055;0.111138317798775,0.419779948680866;0.112144996848305,0.506526488354318;0.109183733379929,0.555451870691800;0.106418083647966,0.607410738054757;0.103261206443177,0.653302635366763;0.0999087155019758,0.696161047653294;0.0925747533317772,0.748412131463658;0.0900603851258762,0.780555940678557;0.0852829302401360,0.825028500960295;0.0807010890908087,0.872534546267508;0.0757280204686558,0.913973621523770;0.0688409535609319,0.949443132211552;0.0578324694774475,0.988423286224159];

%% INTERPOLATE for r/R values
le_R = interp1(tempLE(:,2),tempLE(:,1),r_R,'linear','extrap');
figure(1)
clf(1)
hold on
plot(tempLE(:,2),tempLE(:,1),'-k*')
plot(r_R,le_R,'--dr')
plot(r_R,le_R-c_R,'--db')
ylabel('x-dir (m)')
xlabel('y-dir (m)')
title('Top View of Rotor')
axis equal
grid on
box on
hold off

%% Create input
fcnXMLPANEL(-1*le_R*rotor_radius,r_R*rotor_radius,(zeros(size(r_R,1),1)),c_R*rotor_radius,beta,airfoil_name,output_filename);
