function [pop, max_tip_speed, min_tip_speed] = creation(nvars,num_prop,seed_size,ub,lb,A,b,N_chord,N_dihe,Vars_prop)

totalPopulationSize = seed_size;
pop = zeros(totalPopulationSize,nvars);

A_prop = A(end-3:end,((N_chord + N_dihe)+1):end);
b_prop = b(end-3:end);
lb_prop = lb(((N_chord + N_dihe)+1):end);
ub_prop = ub(((N_chord + N_dihe)+1):end);

air_temp = -1; % Celsius at 8000 ft altitude
c = 331.3*sqrt(1 + air_temp/273.15);
max_tip_speed = 0.84*c; % Mach 0.84 max tip speed
min_tip_speed = 0.5*c;

i = 1;
while i <= totalPopulationSize
    chord = sort([(ub(1:N_chord)-lb(1:N_chord)).*rand(1,N_chord) + lb(1:N_chord)],'descend');
    
    dihedral = [];
    if N_dihe > 0
        dihedral = sort([(ub(N_chord+1:(N_chord + N_dihe))-lb(N_chord+1:(N_chord + N_dihe))).*rand(1,N_dihe) + lb(N_chord+1:(N_chord + N_dihe))],'ascend');
    end
    
    if rand(1) < 0.2
       dihedral = dihedral.*0; 
    end
    
    A_opt = A_prop;
    options = optimoptions('lsqlin','display','none');
    
    prop_rpm = [(ub(N_chord+N_dihe + 2)-lb(N_chord+N_dihe + 2)).*rand(1) + lb(N_chord+N_dihe + 2)];
    
    max_prop_diam = ((max_tip_speed*(60/prop_rpm))/pi)*100;
    if max_prop_diam > ub(N_chord+N_dihe + 1)
        max_prop_diam = ub(N_chord+N_dihe + 1);
    end

    min_prop_diam = ((min_tip_speed*(60/prop_rpm))/pi)*100;
    if min_prop_diam > lb(N_chord+N_dihe + 1)
        min_prop_diam = lb(N_chord+N_dihe + 1);
    end
    
    prop_diam = [(max_prop_diam - min_prop_diam).*rand(1) + min_prop_diam];
    
    lb_prop_temp = lb_prop;
    ub_prop_temp = ub_prop;
    lb_prop_temp(1) = prop_diam;
    ub_prop_temp(1) = prop_diam;
    prop_temp_y = lsqlin(A_opt,b_prop-1,[],[],[],[],lb_prop_temp,ub_prop_temp,[],options); % Solve for remaining chord stations    
    prop_y = prop_temp_y(3:Vars_prop:end)';
        
    pop(i,:) = [chord dihedral prop_diam prop_rpm zeros(1,num_prop*Vars_prop)];
    idx = N_chord+N_dihe + 3;
    for jj = 1:num_prop
        prop_z = [(ub(idx+1)-lb(idx+1)).*rand(1) + lb(idx+1)];
        prop_dir = [(ub(idx+2)-lb(idx+2)).*rand(1) + lb(idx+2)];
        pop(i,idx:idx+Vars_prop-1) = [prop_y(jj) prop_z prop_dir];
        idx = idx + Vars_prop;
    end
    
    if all(A*pop(i,:)' - b <= 0)
        i = i + 1;
    else
       pop(i,:) = pop(i,:).*0; 
    end
end

end

