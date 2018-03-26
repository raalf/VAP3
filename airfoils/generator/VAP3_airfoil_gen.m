function [] = VAP3_airfoil_gen(airfoil_name, coord, Re_range)

ab = dir;
if ~any(cellfun( @(x) strcmp(x,'aux_files'), {ab.name} ))
    mkdir aux_files;
end
    
pol = xfoil(coord, -10:0.25:15, Re_range, 10, [1 1], false);

save(['out_',airfoil_name], 'pol', 'coord', 'Re_range', 'airfoil_name');

end