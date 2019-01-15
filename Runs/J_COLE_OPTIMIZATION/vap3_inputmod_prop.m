function vap3_inputmod_prop(vap_filename, rotor, qmil_output_filename)

str = [];
fp = fopen(vap_filename);
while ~feof(fp)
    str = [str, fgets(fp)];
end
fclose(fp);

flipy = 'FALSE';
if rotor.dir == 0
   flipy = 'TRUE';
   rotor.rpm = -rotor.rpm;
end

str_rotor = sprintf('\n<rotor>\n\t<rpm>%.6f</rpm>\n\t<flipy>%s</flipy>\n\t<collective>%.6f</collective>\n\t<dia>%.6f</dia>\n\t<xhub>%.6f</xhub>\n\t<yhub>%.6f</yhub>\n\t<zhub>%.6f</zhub>\n\t<axisx>%.6f</axisx>\n\t<axisy>%.6f</axisy>\n\t<axisz>%.6f</axisz>\n\t<blades>%d</blades>\n\t<M>%d</M>\n', rotor.rpm, flipy, rotor.collective, rotor.diam, rotor.hub(1), rotor.hub(2), rotor.hub(3), rotor.axis(1), rotor.axis(2), rotor.axis(3), rotor.blades, rotor.m);
rotor_n = 20
str_rotor_geom = fcnQMIL2VAP(qmil_output_filename, rotor_n, 'MH-117');
str_rotor = [str_rotor, str_rotor_geom, sprintf('</rotor>\n')];

k = strfind(str, '</vehicle>') - 6;
str_out = [str(1:k), str_rotor, str(k:end)];

fp = fopen(vap_filename, 'w');
fprintf(fp, str_out);
fclose(fp);

end