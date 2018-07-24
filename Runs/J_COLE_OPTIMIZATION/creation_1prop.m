function pop = creation_1_prop(nvars,seed_size,ub,lb)

totalPopulationSize = seed_size;
pop = zeros(totalPopulationSize,nvars);
for i = 1:totalPopulationSize
    chord = sort([(ub(1:11)-lb(1:11)).*rand(1,11) + lb(1:11)],'descend');
    dihedral = sort([(ub(12:22)-lb(12:22)).*rand(1,11) + lb(12:22)],'ascend');
    prop_diam = [(ub(23)-lb(23)).*rand(1) + lb(23)];
    prop_rpm = [(ub(24)-lb(24)).*rand(1) + lb(24)];
    prop_y = [(ub(25)-lb(25)).*rand(1) + lb(25)];
    prop_z = [(ub(26)-lb(26)).*rand(1) + lb(26)];
    prop_dir = [(ub(27)-lb(27)).*rand(1) + lb(27)];
    pop(i,:) = [chord dihedral prop_diam prop_rpm prop_y prop_z prop_dir];
end

end

