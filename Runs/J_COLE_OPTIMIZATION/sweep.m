clc
clear

addpath('../../')
addpath('../../airfoils')
addpath('./../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../Runs/J_COLE_OPTIMIZATION/')

cores = 37;
parpool(cores,'IdleTimeout',800)
home_dir = pwd;

N_chord = 11;
N_dihe = 0;
N_prop_max = 1;
Vars_prop = 3;

%%
d = linspace(140,200,6);
y = linspace(140, 482*0.95, 12);
z = [-15:5:15];
[D,Y,Z] = meshgrid(d,y,z);
swp_inp = [D(:), Y(:), Z(:)];

J = 77.2/((2250/60)*1.524); % Advance ratio of X57 cruise prop

z = [];
for i = 1:size(swp_inp,1)
   rpm = 60*(77.2/(J*swp_inp(i,1)/100));
   z(i,:) = [linspace(75.6, 53, 11) swp_inp(i,1) rpm swp_inp(i,2:3) 1];
end

%%
max_tip_speed = 277.7821;
min_tip_speed = 165.3465;

parfor i = 1:size(z,1)
    [out, ITER, ITEROUTP] = fcnOBJECTIVE(z(i,:), N_chord, N_dihe, N_prop_max, Vars_prop, max_tip_speed, min_tip_speed, home_dir);
    
    iu_out(i).out = out;
    iu_out(i).ITER = ITER;
    iu_out(i).ITEROUTP = ITEROUTP;
end


