function vap3_inputmod_wing(vap_filename, pan)

str = [];
fp = fopen(vap_filename);
while ~feof(fp)
    str = [str, fgets(fp)];
end
fclose(fp);

wing_M = 1;
str_wing = sprintf(['\n<wing>\n<symmetry>TRUE</symmetry>\n<incidence>0</incidence>\n<trimable>FALSE</trimable>\n<triangular_elements>FALSE</triangular_elements>\n\n<chordwise_elements>', num2str(wing_M), '</chordwise_elements>\n\n<vehicle_x>0</vehicle_x>\n<vehicle_y>0</vehicle_y>\n<vehicle_z>0</vehicle_z>\n']);

str_panel = [];
for ii = 1:size(pan,2)
    if ii == 1
    airfoil = 'FX S 02-196';
    else
    airfoil = 'PSU 94-097';
    end
    
    if ii == 1
        n = 10;
    else
        n = 3;
    end
    str_panel = [str_panel sprintf('\n<panel>\n<spanwise_elements>%d</spanwise_elements>\n<strip_airfoil>%s</strip_airfoil>\n',n,airfoil)];

    for i = 1:size(pan(ii).geom,1)
       str_panel = [str_panel, sprintf('<section>\n\t<wing_x>%.6f</wing_x>\n\t<wing_y>%.6f</wing_y>\n\t<wing_z>%.6f</wing_z>\n\t<chord>%.6f</chord>\n\t<twist>%.6f</twist>\n\t<camber_airfoil>nan</camber_airfoil>\n</section>\n', pan(ii).geom(i,1), pan(ii).geom(i,2), pan(ii).geom(i,3), pan(ii).geom(i,4), pan(ii).geom(i,5))];      
    end

    str_panel = [str_panel, sprintf('</panel>\n\n')];
end



str_wing = [str_wing, str_panel, sprintf('\n</wing>')];

k = strfind(str, '</vehicle>') - 6;
str_out = [str(1:k), str_wing, str(k:end)];

fp = fopen(vap_filename, 'w');
fprintf(fp, str_out);
fclose(fp);

end