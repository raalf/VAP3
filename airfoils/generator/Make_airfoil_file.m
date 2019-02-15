clc
clear

% From a dat file
airfoil_name = 'PSU 94-097';
coord = dlmread([airfoil_name, '.dat'],'',1,0);

VAP3_airfoil_gen(airfoil_name, coord, [150000:200000:8e6], [-2:1:20])

