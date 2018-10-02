function [pop, max_tip_speed, min_tip_speed] = creation_validation(nvars,seed_size,ub,lb,A,b,N_chord,Vars_prop)

totalPopulationSize = seed_size;
pop = zeros(totalPopulationSize,nvars);

air_temp = -1; % Celsius at 8000 ft altitude
c = 331.3*sqrt(1 + air_temp/273.15);
max_tip_speed = 0.84*c; % Mach 0.84 max tip speed
min_tip_speed = 0.5*c;

i = 1;
while i <= totalPopulationSize
    
    idx_seed = round(rand(1,1)*(N_chord-1)) + 1;
    
    A_chord = A(1:end,1:N_chord);
    b_chord = b(1:end);
    lb_chord = lb(1:N_chord);
    ub_chord = ub(1:N_chord); % lower and upper bounds
    lb_chord(idx_seed) = 1; % For the chord distribution problem, as the root chord is set manually
    ub_chord(idx_seed) = 1;
    
    
    % Getting chord distribution, by setting root chord and using least
    % squares to find a chord distribution that does not taper up, and also
    % fits the area constraint
    A_opt = A_chord;
    c_seed = (ub(idx_seed)-lb(idx_seed)).*rand(1,1) + lb(idx_seed); % Randomly set root chord
    A_opt(1,idx_seed) = -c_seed; % We put the known chord in the constriants matrix
    A_opt(end-1,idx_seed) = c_seed;
    A_opt(end,idx_seed) = -c_seed;
    
    options = optimoptions('lsqlin','display','none');
    chord = lsqlin(A_opt,b,[],[],[],[],lb_chord,ub_chord,[],options); % Solve for remaining chord stations
    chord(idx_seed) = c_seed;
    chord = chord(randperm(N_chord));
    prop_diam = [(ub(N_chord + 1)-lb(N_chord + 1)).*rand(1) + lb(N_chord + 1)];

    % C = pi*D
    % C*n = m/s
    % tip speed = pi*D*rps
    % 60*((tip speed)/(pi*D)) = rpm
    
    max_prop_rpm = 60*(max_tip_speed/(pi*prop_diam/100));
    if max_prop_rpm > ub(N_chord + 2)
        max_prop_rpm = ub(N_chord + 2);
    end
    
    min_prop_rpm = 60*(min_tip_speed/(pi*prop_diam/100));
    if min_prop_rpm > ub(N_chord + 2)
        min_prop_rpm = ub(N_chord + 2);
    end
    
    prop_rpm = [(max_prop_rpm - min_prop_rpm).*rand(1) + min_prop_rpm];
    
    prop_y = [(ub(N_chord + 3)-lb(N_chord + 3)).*rand(1) + lb(N_chord + 3)];
    prop_z = [(ub(N_chord + 4)-lb(N_chord + 4)).*rand(1) + lb(N_chord + 4)];
    prop_dir = [(ub(N_chord + 5)-lb(N_chord + 5)).*rand(1) + lb(N_chord + 5)];
    pop(i,:) = [chord' prop_diam prop_rpm prop_y prop_z prop_dir];
    
    if all(A*pop(i,:)' - b <= 0)
        %                 disp(i)
        i = i + 1;
    else
        pop(i,:) = pop(i,:).*0;
    end
end

end

