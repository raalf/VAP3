function pop = creation_1prop(nvars,seed_size,ub,lb)

totalPopulationSize = seed_size;
pop = zeros(totalPopulationSize,nvars);
for i = 1:totalPopulationSize - 2
    chord = sort([(ub(1:11)-lb(1:11)).*rand(1,11) + lb(1:11)],'descend');
    dihedral = sort([(ub(12:22)-lb(12:22)).*rand(1,11) + lb(12:22)],'ascend');
    
    if rand(1) < 0.4
       dihedral = dihedral.*0; 
    end
    
    prop_diam = [(ub(23)-lb(23)).*rand(1) + lb(23)];
    prop_rpm = [(ub(24)-lb(24)).*rand(1) + lb(24)];
    prop_y = [(ub(25)-lb(25)).*rand(1) + lb(25)];
    prop_z = [(ub(26)-lb(26)).*rand(1) + lb(26)];
    prop_dir = [(ub(27)-lb(27)).*rand(1) + lb(27)];
    pop(i,:) = [chord dihedral prop_diam prop_rpm prop_y prop_z prop_dir];
end

pop(totalPopulationSize - 1,:) = [71.3941 68.7533 68.6304 65.9588 64.5861 63.2590 61.8437 61.4030 61.1291 60.0904 59.0279 0 0 0 0 0 0 0 0 0 0 0 160 2226.86 450.979 -12.5485 1];
pop(totalPopulationSize,:) = [85.9495 83.1285 80.1618 77.3493 75.7004 63.8025 60.9942 59.5311 58.3303 57.1046 52.2952 10.054 17.625	17.9299	20.7694	31.5342	33.8535	38.4692	44.6215	49.6259	64.9211	67.2667	160	2097.24	144.366	-13.7057 1];

end

