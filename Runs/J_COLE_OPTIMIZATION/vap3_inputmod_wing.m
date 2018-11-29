function vap3_inputmod_wing(vap_filename, geom)

str = [];
fp = fopen(vap_filename);
while ~feof(fp)
    str = [str, fgets(fp)];
end
fclose(fp);

wing_M = 1
str_wing = sprintf(['\n<wing>\n<incidence>0</incidence>\n<trimable>FALSE</trimable>\n<M>', num2str(wing_M), '</M>\n<xorig>0</xorig>\n<yorig>0</yorig>\n<zorig>0</zorig>\n']);

str_panel = sprintf('\n<panel>\n<airfoil>NASA GA(W)-2 Modified</airfoil>\n<symmetry>1</symmetry>\n<N>2</N>\n');

for i = 1:size(geom,1)
   str_panel = [str_panel, sprintf('<section>\n<x>%.6f</x>\n<y>%.6f</y>\n<z>%.6f</z>\n<chord>%.6f</chord>\n<twist>%.6f</twist>\n</section>\n', geom(i,1), geom(i,2), geom(i,3), geom(i,4), geom(i,5))];      
end

str_panel = [str_panel, sprintf('</panel>\n')];
str_wing = [str_wing, str_panel, sprintf('\n</wing>')];

k = strfind(str, '</vehicle>') - 6;
str_out = [str(1:k), str_wing, str(k:end)];

fp = fopen(vap_filename, 'w');
fprintf(fp, str_out);
fclose(fp);

end