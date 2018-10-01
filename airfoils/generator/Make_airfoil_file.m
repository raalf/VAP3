clc
clear

% From a dat file
airfoil_name = 'NACA 64_2 015';
coord = dlmread([airfoil_name, '.dat'],'',1,0);

VAP3_airfoil_gen(airfoil_name, coord, [150000:150000:8e6], [-10:0.25:20])

