function pop = creation_validation(nvars,seed_size,ub,lb,A,b,N_chord,Vars_prop)

totalPopulationSize = seed_size;
pop = zeros(totalPopulationSize,nvars);

A_chord = A(1:end,1:N_chord);
b_chord = b(1:end);
lb_chord = lb(1:N_chord);
ub_chord = ub(1:N_chord); % lower and upper bounds
lb_chord(1) = 1; % For the chord distribution problem, as the root chord is set manually
ub_chord(1) = 1;
  
i = 1;
while i <= totalPopulationSize
    % Getting chord distribution, by setting root chord and using least
    % squares to find a chord distribution that does not taper up, and also
    % fits the area constraint
    A_opt = A_chord;
    c_root = (ub(1)-lb(1)).*rand(1,1) + lb(1); % Randomly set root chord
    A_opt(1,1) = -c_root; % We put the known chord in the constriants matrix
    A_opt(end-1,1) = c_root;
    A_opt(end,1) = -c_root;
    options = optimoptions('lsqlin','display','none');
    chord = lsqlin(A_opt,b,[],[],[],[],lb_chord,ub_chord,[],options); % Solve for remaining chord stations
    chord(1) = c_root;
    
    dihedral = sort([(ub(12:22)-lb(12:22)).*rand(1,11) + lb(12:22)],'ascend');
    
    prop_diam = [(ub(23)-lb(23)).*rand(1) + lb(23)];
    prop_rpm = [(ub(24)-lb(24)).*rand(1) + lb(24)];
    prop_y = [(ub(25)-lb(25)).*rand(1) + lb(25)];
    prop_z = [(ub(26)-lb(26)).*rand(1) + lb(26)];
    prop_dir = [(ub(27)-lb(27)).*rand(1) + lb(27)];
    pop(i,:) = [chord' dihedral prop_diam prop_rpm prop_y prop_z prop_dir];

    if all(A*pop(i,:)' - b <= 0)
        i = i + 1;
    else
       pop(i,:) = pop(i,:).*0; 
    end
end

end

