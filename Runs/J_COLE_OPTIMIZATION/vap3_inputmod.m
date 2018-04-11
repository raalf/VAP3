function vap_file = vap3_inputmod(rotor)

filename = 'X57_WITHOUT_ROTORS.vap';
str = [];
fp = fopen(filename);
while ~feof(fp)
    str = [str, fgets(fp)];
end
fclose(fp);

str_rotor = sprintf('\n<rotor>\n\t<rpm>%.6f</rpm>\n\t<collective>%.6f</collective>\n\t<dia>%.6f</dia>\n\t<xhub>%.6f</xhub>\n\t<yhub>%.6f</yhub>\n\t<zhub>%.6f</zhub>\n\t<axisx>%.6f</axisx>\n\t<axisy>%.6f</axisy>\n\t<axisz>%.6f</axisz>\n\t<blades>%d</blades>\n\t<M>%d</M>\n', rotor.rpm, rotor.collective, rotor.dia, rotor.hub(1), rotor.hub(2), rotor.hub(3), rotor.axis(1), rotor.axis(2), rotor.axis(3), rotor.blades, rotor.m);
str_rotor_geom = fcnQMIL2VAP('QMIL/out.prop', 10, 'MH-117');
str_rotor = [str_rotor, str_rotor_geom, sprintf('</rotor>\n')];

k = strfind(str, '</wing>') + 7;
str_out = [str(1:k), str_rotor, str(k:end)];

vap_file = [tempname('aux_files'), '.vap'];

fp = fopen(vap_file, 'w');
fprintf(fp, str_out);
fclose(fp);

end