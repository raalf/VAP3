function vap3_inputmod_wing(vap_filename, geom)

str = [];
fp = fopen(vap_filename);
while ~feof(fp)
    str = [str, fgets(fp)];
end
fclose(fp);

str_panel = sprintf('\n<panel>\n<airfoil>NASA GA(W)-2 Modified</airfoil>\n<symmetry>1</symmetry>\n<N>1</N>\n');
% str_panel = '

for i = 1:size(geom,1)
   str_panel = [str_panel, sprintf('<section>\n<x>%.6f</x>\n<y>%.6f</y>\n<z>%.6f</z>\n<chord>%.6f</chord>\n<twist>%.6f</twist>\n</section>\n', geom(i,1), geom(i,2), geom(i,3), geom(i,4), geom(i,5))];      
end

str_panel = [str_panel, sprintf('</panel>\n')];

k = strfind(str, '</wing>') - 6;
str_out = [str(1:k), str_panel, str(k:end)];

fp = fopen(vap_filename, 'w');
fprintf(fp, str_out);
fclose(fp);

end